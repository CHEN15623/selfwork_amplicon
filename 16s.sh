#!/bin/bash
# ==================================================================
# 16S Amplicon Analysis Pipeline  (Linux / conda-adapted)
#
# Original : Windows 11 + R 3.4.5 workflow (USEARCH/VSEARCH/R)
# Adapted  : Ubuntu, single conda environment "16s"
#
#   * ALL steps run inside the conda env "16s" (never in base)
#   * Denoising / taxonomy / alpha-div use the FREE tools that are
#     actually present (the bundled usearch build is the licence-free
#     "personal use" one, which rejects -unoise3, -sintax_summary and
#     -otutab_stats, so vsearch + R replacements are used instead)
#   * Adds FAPROTAX functional prediction with full visualization
#
# Run:  bash /home/cz/16s/16s.sh
# ==================================================================

set -euo pipefail

# === 1. Global configuration =====================================
WD="/home/cz/16s"                 # working directory
DB="${WD}/db"                     # reference databases / scripts
THREADS=8                         # threads for vsearch / mafft
CONDA_ENV="16s"                   # conda environment to use
GROUP_COL="Group"                 # grouping column in metadata
COMPARE_PAIR="M-C"                # differential comparison: case-control
                                  # (metadata groups are C, M, H -> M vs C)
SEED=1                            # random seed
MIN_UNIQUE=10                     # min abundance for a unique/ASV
RAREFY_DEPTH=10000                # even-depth rarefaction

# Display/plot size (mm)
W=89; H=59

# === 2. Activate the conda environment "16s" =====================
# conda is not on PATH by default in a non-login shell, so source the
# hook first. Base environment is deliberately NOT used.
if [ -f /home/cz/miniconda3/etc/profile.d/conda.sh ]; then
    # shellcheck disable=SC1091
    source /home/cz/miniconda3/etc/profile.d/conda.sh
else
    echo "ERROR: conda.sh not found - cannot activate '${CONDA_ENV}'" >&2
    exit 1
fi
conda activate "${CONDA_ENV}"
echo ">>> conda environment : ${CONDA_DEFAULT_ENV}"
echo ">>> R   : $(command -v Rscript)"
echo ">>> python: $(command -v python)"

# keep the licence-free usearch binary available as well
export PATH="${DB}/usearch:$PATH"

# NOTE: every R package used by this pipeline is installed IN the "16s"
# environment itself (amplicon 1.19.0 + agricolae included). No external
# library path is added, so the base conda environment is never touched
# and no R-version mixing can occur.

cd "${WD}"

# === 3. Output directories =======================================
# NOTE: `script/` (relative to the working dir) is required by the bundled
# helper shell scripts (sp_pheatmap.sh / compare_manhattan.sh / ...), which
# write their generated R code into it.
mkdir -p script temp \
    result/raw result/alpha result/beta result/tax \
    result/gg result/compare result/stamp result/lefse \
    result/function result/function/figure

echo ""
echo "=================================================="
echo " 16S pipeline started : $(date '+%F %T')"
echo "=================================================="

echo ">>> [1/11] Environment check ..."
MISSING=""
for t in vsearch mafft FastTree Rscript python seqkit; do
    command -v "$t" >/dev/null 2>&1 || MISSING="${MISSING} ${t}"
done
[ -x "${DB}/usearch/usearch11.0.667_linux32" ] || MISSING="${MISSING} usearch"
if [ -n "${MISSING}" ]; then
    echo "ERROR: missing tools:${MISSING}" >&2
    echo "       install them into the '${CONDA_ENV}' env, e.g.:" >&2
    echo "       conda install -n ${CONDA_ENV} -c conda-forge -c bioconda vsearch mafft fasttree seqkit" >&2
    exit 1
fi
mkdir -p "${DB}/usearch/bin"
ln -sf "${DB}/usearch/usearch11.0.667_linux32" "${DB}/usearch/usearch"
USEARCH="${DB}/usearch/usearch11.0.667_linux32"
echo "    all required tools present"

# Decompress reference databases on demand.
# The repository ships the .gz form to stay within GitHub size limits
# (rdp_16s_v18.fa.gz = 6.5 MB vs 32 MB uncompressed).
for gz in "${DB}/usearch/rdp_16s_v18.fa.gz" "${DB}/gg/97_otus.fa.gz"; do
    out="${gz%.gz}"
    if [ -f "${gz}" ] && [ ! -f "${out}" ]; then
        echo "    decompressing $(basename "${gz}") ..."
        gunzip -c "${gz}" > "${out}"
    fi
done

# Normalise CRLF line endings for every bundled helper script
# (the whole db/script tree was shipped from Windows - CRLF breaks both
#  R and python parsers, so this must run BEFORE anything is executed).
find "${DB}/script" -maxdepth 1 -type f \
    \( -name "*.R" -o -name "*.r" -o -name "*.sh" -o -name "*.pl" -o -name "*.py" \) \
    -exec sed -i 's/\r$//' {} \;
for sub in FAPROTAX_1.2.6 FAPROTAX_1.2.7; do
    [ -d "${DB}/script/${sub}" ] && find "${DB}/script/${sub}" -type f \
        \( -name "*.py" -o -name "*.txt" \) -exec sed -i 's/\r$//' {} \;
done
echo "    CRLF normalisation done"

# === 4. Metadata =================================================
echo ">>> [2/11] Metadata preparation ..."
META="result/metadata.txt"
if [ ! -s "${META}" ]; then
    echo "ERROR: ${META} not found or empty." >&2
    exit 1
fi
sed -i 's/\r$//' "${META}"
# header must be: Sample <tab> Group [<tab> more...]
head -1 "${META}" | grep -q "Sample" || { echo "ERROR: ${META} first column must be 'Sample'" >&2; exit 1; }
SAMPLE_COL=$(head -1 "${META}" | tr '\t' '\n' | grep -n "^Sample$" | cut -d: -f1)
[ -n "${SAMPLE_COL}" ] || { echo "ERROR: no 'Sample' column in ${META}" >&2; exit 1; }
head -1 "${META}" | tr '\t' '\n' | grep -qx "${GROUP_COL}" || { echo "ERROR: no '${GROUP_COL}' column in ${META}" >&2; exit 1; }
tail -n +2 "${META}" | cut -f "${SAMPLE_COL}" > temp/samples.txt
NSAMP=$(wc -l < temp/samples.txt)
echo "    ${NSAMP} samples: $(tr '\n' ' ' < temp/samples.txt)"
GROUP_COL_IDX=$(head -1 "${META}" | tr '\t' '\n' | grep -n "^${GROUP_COL}$" | cut -d: -f1)
echo "    groups: $(tail -n +2 "${META}" | cut -f "${GROUP_COL_IDX}" | sort -u | tr '\n' ' ')"

# === 5. Merge paired-end reads ===================================
echo ">>> [3/11] Merge paired-end reads ..."
rm -f temp/merged_*.fq
EXT=""
for e in fastq.gz fq.gz fastq fq; do
    if [ -f "seq/$(head -1 temp/samples.txt)_1.${e}" ]; then
        EXT="$e"; break
    fi
done
[ -n "${EXT}" ] || { echo "ERROR: cannot find seq/<sample>_1.fastq.gz in seq/" >&2; exit 1; }
echo "    input read extension: .${EXT}"

while read -r s; do
    [ -z "$s" ] && continue
    vsearch --fastq_mergepairs "seq/${s}_1.${EXT}" \
        --reverse "seq/${s}_2.${EXT}" \
        --fastqout "temp/merged_${s}.fq" \
        --fastq_minovlen 20 --fastq_maxdiffs 10 \
        --relabel "${s}." \
        --threads 1 > "temp/merged_${s}.log" 2>&1 &
done < temp/samples.txt
wait
cat temp/merged_*.fq > temp/all.fq
echo "    merged read pairs written to temp/all.fq ($(($(wc -l < temp/all.fq) / 4)) reads)"

# === 6. QC and denoising =========================================
echo ">>> [4/11] Quality control and ASV denoising ..."
# Primers are ALREADY removed in these reads (the V3-V4 amplicon starts
# at the 341F binding site), therefore NO primer trimming is applied.
# Only length / expected-error filtering is performed.
vsearch --fastx_filter temp/all.fq \
    --fastq_minlen 200 --fastq_maxlen 600 \
    --fastq_maxee_rate 0.01 --fastq_maxns 0 \
    --fastaout temp/filtered.fa 2> temp/filter.log
echo "    filtered: $(grep -c '^>' temp/filtered.fa) reads"

vsearch --derep_fulllength temp/filtered.fa \
    --minuniquesize "${MIN_UNIQUE}" --sizeout \
    --relabel Uni_ --output temp/uniques.fa 2> temp/derep.log
echo "    uniques:  $(grep -c '^>' temp/uniques.fa)"

# UNOISE3 denoising - vsearch's free implementation (usearch -unoise3
# is rejected by the licence-free binary).
vsearch --cluster_unoise temp/uniques.fa \
    --minsize "${MIN_UNIQUE}" --unoise_alpha 2.0 \
    --consout temp/zotus.fa 2> temp/unoise.log
# Rename to ASV_N and DROP the "size=" / "seqs=" annotations.
# vsearch --otutabout reports labels WITHOUT any ";size=" suffix, so the
# FASTA must use the very same bare labels or every downstream ID join
# (non-Bacteria filter, -fastx_getseqs, taxonomy) silently matches nothing.
sed 's/^>centroid=/>ASV_/; s/^>zotu/>ASV_/; s/^>Zotu/>ASV_/; s/;.*//' \
    temp/zotus.fa > temp/asv_raw.fa
echo "    ASVs:     $(grep -c '^>' temp/asv_raw.fa)"

# chimera removal against RDP
vsearch --uchime_ref temp/asv_raw.fa \
    --db "${DB}/usearch/rdp_16s_v18.fa" \
    --nonchimeras result/raw/otus.fa \
    --threads "${THREADS}" 2> temp/uchime.log
sed -i 's/\r$//' result/raw/otus.fa
echo "    non-chimeric ASVs: $(grep -c '^>' result/raw/otus.fa)"

# === 7. Feature table + taxonomy =================================
echo ">>> [5/11] Feature table and taxonomy ..."
vsearch --usearch_global temp/filtered.fa \
    --db result/raw/otus.fa --id 0.97 \
    --threads "${THREADS}" --otutabout result/raw/otutab.txt 2> temp/otutab.log
sed -i 's/\r$//' result/raw/otutab.txt
echo "    feature table: $(( $(wc -l < result/raw/otutab.txt) - 1 )) ASVs x $(( $(head -1 result/raw/otutab.txt | wc -w) - 1 )) samples"

vsearch --sintax result/raw/otus.fa \
    --db "${DB}/usearch/rdp_16s_v18.fa" \
    --sintax_cutoff 0.1 --tabbedout result/raw/otus.sintax \
    --threads "${THREADS}" 2> temp/sintax.log
sed -i 's/\r$//' result/raw/otus.sintax

# remove chloroplast / mitochondria / non-bacteria
Rscript "${DB}/script/otutab_filter_nonBac.R" \
    --input result/raw/otutab.txt --taxonomy result/raw/otus.sintax \
    --output result/otutab.txt --stat result/raw/otutab_nonBac.stat \
    --discard result/raw/otus.sintax.discard
echo "    after non-Bacteria filter: $(( $(wc -l < result/otutab.txt) - 1 )) ASVs"

# keep fasta + taxonomy in sync with the filtered table
cut -f 1 result/otutab.txt | tail -n +2 > result/otutab.id
"${USEARCH}" -fastx_getseqs result/raw/otus.fa -labels result/otutab.id \
    -fastaout result/otus.fa > /dev/null 2>&1
awk 'NR==FNR{a[$1]=$0; next} ($1 in a){print a[$1]}' \
    result/raw/otus.sintax result/otutab.id > result/otus.sintax

# === 8. Rarefaction ==============================================
echo ">>> [6/11] Even-depth rarefaction (depth=${RAREFY_DEPTH}) ..."
Rscript "${DB}/script/otutab_rare.R" --input result/otutab.txt \
    --depth "${RAREFY_DEPTH}" --seed "${SEED}" \
    --normalize result/otutab_rare.txt --output result/alpha/vegan.txt

# === 9. Alpha diversity ==========================================
echo ">>> [7/11] Alpha diversity ..."
"${USEARCH}" -alpha_div result/otutab_rare.txt -output result/alpha/alpha.txt > /dev/null 2>&1
sed -i '1s/^Sample/Sample/' result/alpha/alpha.txt

# rarefaction curve data (vegan based, usearch -alpha_div_rare is paid)
Rscript "${DB}/script/alpha_rare_curve_gen.R" \
    -i result/otutab_rare.txt -o result/alpha/alpha_rare.txt

# boxplot per index
for idx in $(head -n1 result/alpha/vegan.txt | cut -f 2-); do
    Rscript "${DB}/script/alpha_boxplot.R" --alpha_index "${idx}" \
        --input result/alpha/vegan.txt --design "${META}" \
        --group "${GROUP_COL}" --output result/alpha/ --width "${W}" --height "${H}"
done
# barplot
# NOTE: amplicon::alpha_barplot() calls the defunct dplyr verb
# `summarise_()` and always errors with dplyr >= 1.1, so the maintained
# drop-in replacement alpha_barplot_fix.R is used instead (identical CLI).
Rscript "${DB}/script/alpha_barplot_fix.R" --alpha_index richness \
    --input result/alpha/vegan.txt --design "${META}" \
    --group "${GROUP_COL}" --output result/alpha/ --width "${W}" --height "${H}"

# rarefaction curve figure
Rscript "${DB}/script/alpha_rare_curve_custom.R" \
    --input result/otutab_rare.txt --design "${META}" \
    --group "${GROUP_COL}" --output result/alpha/ --width 120 --height "${H}"

# group-mean table + Venn
Rscript "${DB}/script/otu_mean.R" --input result/otutab.txt --metadata "${META}" \
    --group "${GROUP_COL}" --thre 0 --scale TRUE --zoom 100 --all TRUE \
    --type mean --output result/otutab_mean.txt
ALLGROUPS=()
while IFS= read -r g; do [ -n "$g" ] && ALLGROUPS+=("$g"); done < <(tail -n +2 "${META}" | cut -f 2 | sort -u)
awk 'BEGIN{OFS=FS="\t"}
     FNR==1{for(i=2;i<=NF;i++) a[i]=$i; print "OTU","Group"; next}
     {for(i=2;i<=NF;i++) if($i+0>0.05) print $1, a[i]}' \
    result/otutab_mean.txt > result/alpha/otu_group_exist.txt
if [ "${#ALLGROUPS[@]}" -ge 2 ] && [ -f "${DB}/script/sp_vennDiagram.sh" ]; then
    VARG=""; NG=0; VNAMES=()
    for g in "${ALLGROUPS[@]}"; do
        [ "${NG}" -ge 4 ] && break
        # skip the combined "All" column produced by otu_mean.R
        [ "${g}" == "All" ] && continue
        NG=$((NG+1))
        VNAMES+=("${g}")
        letter=$(printf "\\$(printf '%03o' $((96+NG)))")
        VARG="${VARG} -${letter} ${g}"
    done
    if [ "${NG}" -ge 2 ]; then
        VPNAME=$(printf '%s_' "${VNAMES[@]}")
        VPNAME=${VPNAME%_}
        echo "    Venn groups: ${VNAMES[*]}"
        bash "${DB}/script/sp_vennDiagram.sh" -f result/alpha/otu_group_exist.txt \
            ${VARG} -w 3 -u 3 -p "${VPNAME}" || echo "    [warn] Venn failed (skipped)"
    else
        echo "    [skip] Venn needs >=2 groups"
    fi
else
    echo "    [skip] Venn needs >=2 groups"
fi

# === 10. Beta diversity ==========================================
echo ">>> [8/11] Beta diversity ..."
mafft --thread "${THREADS}" --quiet result/otus.fa > temp/otus_aligned.fa 2>/dev/null
FastTree -nt -quiet temp/otus_aligned.fa > result/otus.tree 2>/dev/null
echo "    tree built: $(grep -c '>' result/otus.fa) tips"

Rscript "${DB}/script/beta_div_gen.R" \
    --input result/otutab_rare.txt --tree result/otus.tree --prefix result/beta/

for label_flag in FALSE TRUE; do
    out_suffix=$( [ "$label_flag" == "TRUE" ] && echo ".label" || echo "" )
    Rscript "${DB}/script/beta_pcoa.R" --input result/beta/bray_curtis.txt \
        --design "${META}" --group "${GROUP_COL}" --label "${label_flag}" \
        --width "${W}" --height "${H}" \
        --output "result/beta/bray_curtis.pcoa${out_suffix}.pdf"
done

Rscript "${DB}/script/beta_cpcoa.R" --input result/beta/bray_curtis.txt \
    --design "${META}" --group "${GROUP_COL}" --label TRUE \
    --width "${W}" --height "${H}" --output result/beta/bray_curtis.cpcoa.label.pdf

cut -f 1-2 "${META}" > temp/group.txt
if [ -f "${DB}/script/sp_pheatmap.sh" ]; then
    bash "${DB}/script/sp_pheatmap.sh" -f result/beta/bray_curtis.txt -H 'TRUE' \
        -u 6.9 -v 5.6 -P temp/group.txt -Q temp/group.txt \
        || echo "    [warn] beta heatmap failed (skipped)"
fi

# === 11. Taxonomy composition ====================================
echo ">>> [9/11] Taxonomy summary and figures ..."
# 8-field taxonomy table
cut -f 1,4 result/otus.sintax | sed 's/\td/\tk/;s/:/__/g;s/,/;/g;s/"//g' > result/taxonomy2.txt
awk 'BEGIN{OFS=FS="\t"}
     {delete a; a["k"]="Unassigned";a["p"]="Unassigned";a["c"]="Unassigned";
      a["o"]="Unassigned";a["f"]="Unassigned";a["g"]="Unassigned";a["s"]="Unassigned";
      split($2,x,";"); for(i in x){split(x[i],b,"__"); a[b[1]]=b[2];}
      print $1,a["k"],a["p"],a["c"],a["o"],a["f"],a["g"],a["s"];}' \
    result/taxonomy2.txt > temp/otus.tax
sed 's/;/\t/g; s/.__//g;' temp/otus.tax | cut -f 1-8 | \
    sed '1 s/^/OTUID\tKingdom\tPhylum\tClass\tOrder\tFamily\tGenus\tSpecies\n/' \
    > result/taxonomy.txt

# rank-level relative-abundance tables (free replacement for
# `usearch -sintax_summary`, which the licence-free build rejects)
Rscript "${DB}/script/sintax_rank_summary.R" \
    --input result/otutab_rare.txt --taxonomy result/otus.sintax \
    --output result/tax/ --ranks p,c,o,f,g

for rank in p c o f g; do
    [ -s "result/tax/sum_${rank}.txt" ] || { echo "    [skip] sum_${rank}.txt missing"; continue; }
    Rscript "${DB}/script/tax_stackplot.R" --input "result/tax/sum_${rank}.txt" \
        --design "${META}" --group "${GROUP_COL}" --legend 8 \
        --output "result/tax/sum_${rank}.stackplot" --width "${W}" --height "${H}"

    if [ "${rank}" == "c" ]; then
        Rscript "${DB}/script/tax_circlize.R" --input "result/tax/sum_${rank}.txt" \
            --design "${META}" --group "${GROUP_COL}" --legend 5 \
            || echo "    [warn] circlize failed"
        mv -f circlize.pdf "result/tax/sum_${rank}.circlize.pdf" 2>/dev/null || true
        mv -f circlize_legend.pdf "result/tax/sum_${rank}.circlize_legend.pdf" 2>/dev/null || true
    fi
done

# Treemap of the most abundant ASVs.
# amplicon::format2maptree() needs the heavy Bioconductor stack
# (phyloseq + ggtree + ggraph). If it is unavailable we fall back to the
# self-contained treemap implementation, which needs only the `treemap`
# package, so this step always produces a figure.
if Rscript -e 'quit(status = !requireNamespace("phyloseq", quietly=TRUE))' 2>/dev/null; then
    Rscript "${DB}/script/tax_maptree.R" --input result/otutab.txt \
        --taxonomy result/taxonomy.txt --output result/tax/tax_maptree.pdf \
        --topN 100 --width 183 --height 118 \
        || Rscript "${DB}/script/tax_maptree_treemap.R" --input result/otutab.txt \
             --taxonomy result/taxonomy.txt --output result/tax/tax_maptree.pdf \
             --topN 100 --width 183 --height 118
else
    echo "    phyloseq not installed - using treemap fallback"
    Rscript "${DB}/script/tax_maptree_treemap.R" --input result/otutab.txt \
        --taxonomy result/taxonomy.txt --output result/tax/tax_maptree.pdf \
        --topN 100 --width 183 --height 118
fi

# === 12. GreenGenes reference mapping ============================
echo ">>> [10/11] GreenGenes 97% reference mapping ..."
vsearch --usearch_global temp/filtered.fa --db "${DB}/gg/97_otus.fa" \
    --id 0.97 --threads "${THREADS}" --otutabout result/gg/otutab.txt 2> temp/gg.log
{
    echo "GreenGenes 97% OTU table"
    echo "OTUs   : $(( $(wc -l < result/gg/otutab.txt) - 1 ))"
    echo "Samples: $(( $(head -1 result/gg/otutab.txt | wc -w) - 1 ))"
} > result/gg/otutab.stat
cat result/gg/otutab.stat

# === 13. Functional prediction (FAPROTAX) ========================
echo ">>> [11/11] FAPROTAX functional prediction ..."
FAP_DIR="${DB}/script/FAPROTAX_1.2.7"
mkdir -p result/function/figure
if [ -f "${FAP_DIR}/collapse_table.py" ] && [ -f "${FAP_DIR}/FAPROTAX.txt" ]; then
    # The bundled 1.2.7 script is not compatible with numpy>=2
    # (`from numpy import NaN` was removed) - build a patched copy in
    # temp/ so the original database stays untouched.
    cp -f "${FAP_DIR}/collapse_table.py" temp/collapse_table.py
    sed -i 's/\r$//' temp/collapse_table.py
    sed -i 's/^from numpy import NaN$/from numpy import nan as NaN/' temp/collapse_table.py

    # FAPROTAX input: classic table with row names in an "OTUID" column and
    # taxonomic lineage in a "taxonomy" column (see FAPROTAX README).
    awk 'BEGIN{OFS="\t"}
         NR==FNR{n[$1]=$4; next}
         FNR==1{print "OTUID", substr($0, index($0,"\t")+1), "taxonomy"; next}
         {print $1, substr($0, index($0,"\t")+1), (($1 in n) ? n[$1] : "Unassigned")}' \
        result/otus.sintax result/otutab.txt > temp/fap_otu.txt
    echo "    FAPROTAX input: $(( $(wc -l < temp/fap_otu.txt) - 1 )) ASVs"

    python temp/collapse_table.py \
        -i temp/fap_otu.txt -g "${FAP_DIR}/FAPROTAX.txt" \
        -o result/function/faprotax.txt \
        -r result/function/faprotax_report.txt \
        -d taxonomy --omit_columns 0 \
        --column_names_are_in first_data_line \
        --normalize_collapsed columns_before_collapsing \
        2> result/function/faprotax.log || {
            echo "    [warn] FAPROTAX failed - see result/function/faprotax.log"
        }

    if [ -s result/function/faprotax.txt ]; then
        echo "    functional groups: $(( $(wc -l < result/function/faprotax.txt) - 1 ))"
        # visualization
        Rscript "${DB}/script/function_faprotax_plot.R" \
            --input result/function/faprotax.txt --design "${META}" \
            --group "${GROUP_COL}" --compare "${COMPARE_PAIR}" \
            --output result/function/ --legend 12 --topN 30 \
            --width 120 --height 90 || echo "    [warn] FAPROTAX plotting failed"
    else
        echo "    [warn] FAPROTAX produced no output"
    fi
else
    echo "    [skip] FAPROTAX database not found under ${FAP_DIR}"
fi

# === 14. Differential analysis / STAMP / LEfSe ===================
echo ">>> [12/12] Differential analysis (compare=${COMPARE_PAIR}) ..."
set +e
Rscript "${DB}/script/compare.R" --input result/otutab.txt --design "${META}" \
    --group "${GROUP_COL}" --compare "${COMPARE_PAIR}" \
    --threshold 0.1 --method edgeR --output result/compare/
RC=$?
set -e
if [ ${RC} -eq 0 ] && [ -s "result/compare/${COMPARE_PAIR}.txt" ]; then
    Rscript "${DB}/script/compare_volcano.R" \
        --input "result/compare/${COMPARE_PAIR}.txt" \
        --output "result/compare/${COMPARE_PAIR}.volcano.pdf" \
        --width "${W}" --height "${H}" || echo "    [warn] volcano failed"

    # Detect how many leading columns of the comparison table hold statistics
    # rather than samples. compare_heatmap.sh needs this offset via -l, and a
    # wrong value makes it select no sample columns at all (pheatmap then
    # aborts with "'gpar' component 'fill' has length zero").
    # The table only contains the samples of the two compared groups, so the
    # offset is found by locating the first header field that is a sample ID.
    NON_DATA_LINE=$(head -1 "result/compare/${COMPARE_PAIR}.txt" | \
        awk -v sfile=temp/samples.txt \
            'BEGIN{while((getline l < sfile)>0){s[l]=1; s["X" l]=1}}
             {for(i=1;i<=NF;i++) if($i in s){print i-1; exit}}')
    [ -n "${NON_DATA_LINE}" ] || NON_DATA_LINE=7
    echo "    compare heatmap: ${NON_DATA_LINE} annotation column(s)"
    [ -f "${DB}/script/compare_heatmap.sh" ] && \
        bash "${DB}/script/compare_heatmap.sh" -i "result/compare/${COMPARE_PAIR}.txt" \
            -l "${NON_DATA_LINE}" \
            -d "${META}" -A "${GROUP_COL}" -t result/taxonomy.txt -w 8 -h 5 -s 7 \
            -o "result/compare/${COMPARE_PAIR}" || echo "    [warn] compare heatmap skipped"

    if [ -f "${DB}/script/compare_manhattan.sh" ]; then
        bash "${DB}/script/compare_manhattan.sh" -i "result/compare/${COMPARE_PAIR}.txt" \
            -t result/taxonomy.txt -p result/tax/sum_p.txt -w 250 -v "${H}" -s 7 -l 10 \
            -o "result/compare/${COMPARE_PAIR}.manhattan.p.pdf" || true
        bash "${DB}/script/compare_manhattan.sh" -i "result/compare/${COMPARE_PAIR}.txt" \
            -t result/taxonomy.txt -p result/tax/sum_g.txt -w 250 -v "${H}" -s 7 -l 10 \
            -L Genus -o "result/compare/${COMPARE_PAIR}.manhattan.g.pdf" || true
    fi
else
    echo "    [warn] compare.R failed - differential plots skipped"
fi

# STAMP
Rscript "${DB}/script/format2stamp.R" --input result/otutab.txt \
    --taxonomy result/taxonomy.txt --threshold 0.01 --output result/stamp/tax \
    || echo "    [warn] format2stamp skipped"
# NOTE: --pvalue is the significance cut-off and --threshold the minimum
# mean relative abundance. With only 3 replicates per group nothing
# usually survives FDR, in which case compare_stamp.R exits cleanly and
# simply writes no figure.
if [ -s result/stamp/tax_6Genus.txt ]; then
    Rscript "${DB}/script/compare_stamp.R" --input result/stamp/tax_6Genus.txt \
        --metadata "${META}" --group "${GROUP_COL}" --compare "${COMPARE_PAIR}" \
        --threshold 0.001 --pvalue 0.05 --method "t.test" --width 280 --height 159 \
        --output "result/stamp/${COMPARE_PAIR}" || echo "    [warn] STAMP plot skipped"
fi

# LEfSe
Rscript "${DB}/script/format2lefse.R" --input result/otutab.txt \
    --taxonomy result/taxonomy.txt --design "${META}" \
    --group "${GROUP_COL}" --threshold 0.4 --output result/lefse/LEfSe \
    || echo "    [warn] format2lefse skipped"

# === 15. Clean-up and checksums ==================================
echo ">>> clean-up ..."
find temp -maxdepth 1 -type f \( -name "*.fq" -o -name "*.fa" \) -delete 2>/dev/null || true

( cd seq && md5sum *_1."${EXT}" > md5sum1.txt && md5sum *_2."${EXT}" > md5sum2.txt && \
  paste md5sum1.txt md5sum2.txt | awk '{print $2"\t"$1"\t"$4"\t"$3}' | sed 's/*//g' > ../result/md5sum.txt && \
  rm -f md5sum1.txt md5sum2.txt ) || echo "    [warn] md5 step skipped"

echo ""
echo "=================================================="
echo "  16S analysis finished : $(date '+%F %T')"
echo "  Key results:"
echo "    result/otutab.txt            feature table"
echo "    result/taxonomy.txt          taxonomy"
echo "    result/alpha/                alpha diversity"
echo "    result/beta/                 beta diversity / PCoA"
echo "    result/tax/                  composition figures"
echo "    result/function/             FAPROTAX prediction"
echo "    result/compare/              differential analysis"
echo "=================================================="
