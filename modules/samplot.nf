process SAMPLOT {
    tag "$meta.id"
    label 'process_low'
    
    conda 'bioconda::samplot=1.3.0'

    input:
    tuple val(meta), path(vcf), path(bam), path(bai)
    path reference
    val max_plots
    val caller
    
    output:
    tuple val(meta), path("*.png"), emit: plots, optional: true
    path "*_mqc.json", emit: multiqc, optional: true
    
    script:
    def prefix = "${meta.id}_${caller}"
    """
    set -euo pipefail

    # Parse VCF and write regions file
    python3 << 'PYCODE'
import gzip
import json
import base64

def smart_open(filename):
    with open(filename, 'rb') as f:
        magic = f.read(2)
        f.seek(0)
        if magic == b'\\x1f\\x8b':
            return gzip.open(filename, 'rt')
        else:
            return open(filename, 'r')

VCF_PATH   = "${vcf}"
MAX_PLOTS  = ${max_plots}
REGIONS_FN = "samplot_regions.tsv"

regions = []
line_num = 0
skipped = 0

with smart_open(VCF_PATH) as f:
    for line in f:
        if line.startswith("#"):
            continue
        
        line_num += 1
        
        if len(regions) >= MAX_PLOTS:
            break

        fields = line.strip().split("\\t")
        if len(fields) < 8:
            print(f"SKIPPED line {line_num}: not enough fields ({len(fields)})")
            skipped += 1
            continue

        chrom = fields[0]
        try:
            pos = int(fields[1])
        except ValueError:
            print(f"SKIPPED line {line_num}: invalid position {fields[1]}")
            skipped += 1
            continue

        info = dict(
            item.split("=", 1) if "=" in item else (item, True)
            for item in fields[7].split(";") if item
        )

        svtype = info.get("SVTYPE", "UNK")

        svlen_raw = info.get("SVLEN")
        if svlen_raw is None:
            svlen_val = 1000
        else:
            try:
                svlen_val = abs(int(svlen_raw))
            except ValueError:
                print(f"SKIPPED line {line_num}: invalid SVLEN {svlen_raw}")
                skipped += 1
                svlen_val = 1000

        # Only plot interesting SVs (>500bp)
        if svlen_val < 500:
            print(f"SKIPPED line {line_num}: SVLEN {svlen_val} < 500bp")
            skipped += 1
            continue

        start = max(0, pos - 500)
        end   = pos + svlen_val + 500

        regions.append({
            'chrom': chrom,
            'start': start,
            'end': end,
            'pos': pos,
            'svtype': svtype
        })
        print(f"ADDED region {len(regions)}: {chrom}:{pos} ({svtype})")

print(f"\\nTotal lines processed: {line_num}")
print(f"Total skipped: {skipped}")
print(f"Total regions added: {len(regions)}")

# Write regions file for bash processing
with open(REGIONS_FN, "w") as out:
    for r in regions:
        out.write(f"{r['chrom']}\\t{r['start']}\\t{r['end']}\\t{r['pos']}\\t{r['svtype']}\\n")

# Save regions for later MultiQC processing
with open('regions.json', 'w') as f:
    json.dump(regions, f)
PYCODE

    # If no regions, exit
    if [ ! -s samplot_regions.tsv ]; then
        echo "No regions to plot"
        exit 0
    fi

    # Generate plots
    plot_num=0
    : > plotted_regions.tsv
    while IFS=\$'\\t' read -r chrom start end pos svtype; do
        [ -z "\$chrom" ] && continue

        out="${prefix}.\${plot_num}_\${chrom}_\${pos}_\${svtype}.png"
        echo "Plotting \$chrom:\$pos (\$svtype) -> \$out"

        if samplot plot \\
            -c "\$chrom" \\
            -s "\$start" \\
            -e "\$end" \\
            -b "${bam}" \\
            -o "\$out" \\
            -t "\$svtype" \\
            -r "${reference}"; then
            echo "SUCCESS: Plot \$((plot_num + 1)) created"
            printf '%s\\t%s\\t%s\\t%s\\n' "\$out" "\$chrom" "\$pos" "\$svtype" >> plotted_regions.tsv
            plot_num=\$((plot_num + 1))
        else
            echo "FAILED: Could not create plot for \$chrom:\$pos (\$svtype)"
        fi
    done < samplot_regions.tsv
    echo "Total plots successfully created: \$plot_num"

    # Create MultiQC custom content with EMBEDDED images
    python3 << 'PYCODE2'
import json
import base64

# Only plots that were created successfully are listed here, so each picture keeps its own label
plotted = []
with open('plotted_regions.tsv') as f:
    for row in f:
        png_file, chrom, pos, svtype = row.rstrip("\\n").split("\\t")
        plotted.append((png_file, {'chrom': chrom, 'pos': pos, 'svtype': svtype}))

# Create HTML content with embedded images
html_content = '<div style="max-width: 100%;">'
html_content += '<h4>Samplot Visualizations</h4>'

for png_file, region in plotted:
    with open(png_file, 'rb') as f:
        img_data = base64.b64encode(f.read()).decode('utf-8')
    
    html_content += f'''
    <div style="margin-bottom: 20px; border: 1px solid #ddd; padding: 10px;">
        <h5>{region['svtype']} at {region['chrom']}:{region['pos']}</h5>
        <img src="data:image/png;base64,{img_data}" style="width: 100%; max-width: 1200px;" />
    </div>
    '''

html_content += '</div>'

mqc_data = {
    "id": "samplot_${caller}",
    "section_name": "Samplot Visualizations - ${caller}",
    "description": f"Visual inspection of {len(plotted)} structural variants for ${meta.id}",
    "plot_type": "html",
    "data": html_content
}

with open('${prefix}_samplot_mqc.json', 'w') as f:
    json.dump(mqc_data, f, indent=2)
PYCODE2
    """
}