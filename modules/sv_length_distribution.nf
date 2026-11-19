process SV_LENGTH_DISTRIBUTION {
    tag "${meta.id}_${caller}"
    label 'process_single'
    
    conda 'python>=3.9 matplotlib>=3.5 seaborn>=0.12 pandas>=1.3 numpy>=1.20 pyyaml plotly'
    
    input:
    tuple val(meta), path(vcf)
    val caller
    
    output:
    tuple val(meta), path("*.length_dist.txt"), emit: stats
    path "*_mqc.json", emit: multiqc
    
    script:
    def prefix = "${meta.id}_${caller}"
    """
    #!/usr/bin/env python3
    import gzip
    import json
    
    # Helper function to open both gzipped and uncompressed files
    def smart_open(filename):
        with open(filename, 'rb') as f:
            magic = f.read(2)
            f.seek(0)
            if magic == b'\\x1f\\x8b':
                return gzip.open(filename, 'rt')
            else:
                return open(filename, 'r')
    
    # Parse VCF for SV lengths
    sv_types = {'DEL': [], 'INS': [], 'DUP': [], 'INV': [], 'BND': []}
    
    with smart_open('${vcf}') as f:
        for line in f:
            if line.startswith('#'):
                continue
            fields = line.strip().split('\\t')
            if len(fields) < 8:
                continue
                
            info = dict(item.split('=', 1) if '=' in item else (item, True) 
                       for item in fields[7].split(';') if item)
            
            svtype = info.get('SVTYPE', 'UNKNOWN')
            svlen_raw = info.get('SVLEN', '0')
            
            try:
                svlen = abs(int(svlen_raw))
            except (ValueError, TypeError):
                svlen = 0
            
            if svtype in sv_types and svlen > 0:
                sv_types[svtype].append(svlen)
    
    # Write stats
    with open('${prefix}.length_dist.txt', 'w') as out:
        out.write('SV_Type\\tCount\\tMin\\tMax\\tMean\\tMedian\\n')
        for svtype, lengths in sv_types.items():
            if lengths:
                out.write(f'{svtype}\\t{len(lengths)}\\t{min(lengths)}\\t{max(lengths)}\\t'
                         f'{sum(lengths)/len(lengths):.1f}\\t{sorted(lengths)[len(lengths)//2]}\\n')
    
    # Create MultiQC custom content with INTERACTIVE BAR PLOT
    # Prepare data for bargraph
    plot_data = {}
    categories = {}
    
    for svtype, lengths in sv_types.items():
        if lengths:
            plot_data['${prefix}'] = plot_data.get('${prefix}', {})
            plot_data['${prefix}'][svtype] = len(lengths)
            categories[svtype] = {'name': svtype}
    
    mqc_data = {
        'id': 'sv_length_dist_${caller}',
        'section_name': 'SV Counts - ${caller}',
        'description': 'Number of structural variants by type for ${meta.id}',
        'plot_type': 'bargraph',
        'pconfig': {
            'id': 'sv_count_bargraph_${caller}',
            'title': 'SV Counts by Type - ${caller}',
            'ylab': 'Count',
            'cpswitch_counts_label': 'Number of SVs'
        },
        'data': plot_data
    }
    
    # Also add a table with statistics
    table_data = {}
    for svtype, lengths in sv_types.items():
        if lengths:
            table_data[svtype] = {
                'Count': len(lengths),
                'Min (bp)': min(lengths),
                'Max (bp)': max(lengths),
                'Mean (bp)': round(sum(lengths) / len(lengths), 1),
                'Median (bp)': sorted(lengths)[len(lengths)//2]
            }
    
    if table_data:
        mqc_table = {
            'id': 'sv_stats_table_${caller}',
            'section_name': 'SV Statistics - ${caller}',
            'description': 'Detailed statistics for structural variants in ${meta.id}',
            'plot_type': 'table',
            'pconfig': {
                'id': 'sv_stats_table_${caller}',
                'title': 'SV Length Statistics - ${caller}'
            },
            'data': {'${prefix}': table_data}
        }
        
        # Combine both sections
        combined = {
            'id': 'sv_analysis_${caller}',
            'section_name': 'SV Analysis - ${caller}',
            'plot_type': 'bargraph',
            'data': plot_data,
            'pconfig': mqc_data['pconfig']
        }
    
    with open('${prefix}_sv_mqc.json', 'w') as f:
        json.dump(mqc_data, f, indent=2)
    """
}