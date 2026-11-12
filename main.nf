#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

// Import modules
include { MINIMAP2   } from './modules/minimap2'
include { PEPPER     } from './modules/pepper'
include { CUTESV     } from './modules/cutesv'
include { CLAIR3     } from './modules/clair3'
include { SNIFFLES2  } from './modules/sniffles2'
include { MULTIQC    } from './modules/multiqc'
include { BCFTOOLS_STATS as BCFTOOLS_STATS_PEPPER   } from './modules/bcftools_stats'
include { BCFTOOLS_STATS as BCFTOOLS_STATS_CUTESV   } from './modules/bcftools_stats'
include { BCFTOOLS_STATS as BCFTOOLS_STATS_CLAIR3   } from './modules/bcftools_stats'
include { BCFTOOLS_STATS as BCFTOOLS_STATS_SNIFFLES } from './modules/bcftools_stats'
include { RTG_VCFSTATS as RTG_VCFSTATS_PEPPER   } from './modules/rtg_vcfstats'
include { RTG_VCFSTATS as RTG_VCFSTATS_CUTESV   } from './modules/rtg_vcfstats'
include { RTG_VCFSTATS as RTG_VCFSTATS_CLAIR3   } from './modules/rtg_vcfstats'
include { RTG_VCFSTATS as RTG_VCFSTATS_SNIFFLES } from './modules/rtg_vcfstats'
include { SV_LENGTH_DISTRIBUTION as SV_LENGTH_DISTRIBUTION_CUTESV } from './modules/sv_length_distribution'
include { SV_LENGTH_DISTRIBUTION as SV_LENGTH_DISTRIBUTION_SNIFFLES } from './modules/sv_length_distribution'
include { SAMPLOT as SAMPLOT_CUTESV } from './modules/samplot'
include { SAMPLOT as SAMPLOT_SNIFFLES } from './modules/samplot'

def validateParams() {
    def errors = []
    
    if (!params.input) errors << "--input is required (path to fastq files)"
    if (!params.outdir) errors << "--outdir is required (output directory)"
    if (!params.reference) errors << "--reference is required"
    if (!params.minimap2_preset) errors << "--minimap2_preset is required (e.g., 'lr:hq')"
    if (params.run_clair3 && !params.clair3_model) errors << "--clair3_model required when --run_clair3=true"
    if (params.run_clair3 && !params.clair3_platform) errors << "--clair3_platform required when --run_clair3=true (e.g., 'ont', 'hifi')"
    if (params.run_pepper && !params.pepper_preset) errors << "--pepper_preset required when --run_pepper=true (e.g., '--ont_r10_q20')"
    
    if (errors) error "Validation failed:\n  - ${errors.join('\n  - ')}"
}

workflow {
    validateParams()    
    // 1. Create input channel with metadata
    // Meta map contains sample information that travels with the data
    Channel
        .fromPath(params.input, checkIfExists: true)
        .map { fastq ->
            def meta = [
                id: fastq.baseName,  // Sample ID from filename
            ]
            [ meta, fastq ]
        }
        .set { ch_input }
    
    // 2. Prepare reference files
    // Single files shared by every sample MUST be value channels. A queue
    // channel (Channel.fromPath) holds one item, is consumed by the first
    // sample, and every other sample is silently dropped.
    ch_reference     = Channel.value(file(params.reference, checkIfExists: true))
    ch_reference_fai = Channel.value(file(params.reference + '.fai', checkIfExists: true))

    // 3. Handle optional inputs
    // Tandem repeats for Sniffles2 (OPTIONAL)
    ch_tandem_repeats = params.tandem_repeats
        ? Channel.value(file(params.tandem_repeats, checkIfExists: true))
        : Channel.value([])  // Empty value if not provided

    // Clair3 model directory (required if running Clair3)
    ch_clair3_model = params.run_clair3
        ? Channel.value(file(params.clair3_model, checkIfExists: true))
        : Channel.value([])
    
    // 4. Run alignment with Minimap2
    MINIMAP2(
        ch_input,
        ch_reference
    )
    
    // 5. Run variant callers (conditional based on params)
    
    // PEPPER - runs only if params.run_pepper == true
    PEPPER(
        MINIMAP2.out.bam.join(MINIMAP2.out.bai),
        ch_reference,
        ch_reference_fai
    )
    
    // cuteSV - runs only if params.run_cutesv == true
    CUTESV(
        MINIMAP2.out.bam.join(MINIMAP2.out.bai),
        ch_reference
    )
    
    // Clair3 - runs only if params.run_clair3 == true
    CLAIR3(
        MINIMAP2.out.bam.join(MINIMAP2.out.bai),
        ch_reference,
        ch_reference_fai,
        ch_clair3_model
    )
    
    // Sniffles2 - runs only if params.run_sniffles2 == true
    // Handles optional tandem_repeats input
    SNIFFLES2(
        MINIMAP2.out.bam.join(MINIMAP2.out.bai),
        ch_reference,
        ch_reference_fai,
        ch_tandem_repeats
    )
    
    // 6. Stats/QC 
    ch_pepper_stats             = BCFTOOLS_STATS_PEPPER(PEPPER.out.vcf, 'pepper')
    ch_pepper_rtg               = RTG_VCFSTATS_PEPPER(PEPPER.out.vcf, 'pepper')

    ch_cutesv_stats             = BCFTOOLS_STATS_CUTESV(CUTESV.out.vcf, 'cutesv')
    ch_cutesv_rtg               = RTG_VCFSTATS_CUTESV(CUTESV.out.vcf, 'cutesv')
    ch_cutesv_length_dist       = SV_LENGTH_DISTRIBUTION_CUTESV(CUTESV.out.vcf, 'cutesv')
    ch_cutesv_samplots          = SAMPLOT_CUTESV(CUTESV.out.vcf.join(MINIMAP2.out.bam).join(MINIMAP2.out.bai), ch_reference, params.samplot_max_plots, 'cutesv')

    ch_clair3_stats             = BCFTOOLS_STATS_CLAIR3(CLAIR3.out.vcf, 'clair3')
    ch_clair3_rtg               = RTG_VCFSTATS_CLAIR3(CLAIR3.out.vcf, 'clair3')

    ch_sniffles2_stats          = BCFTOOLS_STATS_SNIFFLES(SNIFFLES2.out.vcf, 'sniffles2')
    ch_sniffles2_rtg            = RTG_VCFSTATS_SNIFFLES(SNIFFLES2.out.vcf, 'sniffles2')
    ch_sniffles2_length_dist    = SV_LENGTH_DISTRIBUTION_SNIFFLES(SNIFFLES2.out.vcf, 'sniffles2')
    ch_sniffles2_samplots       = SAMPLOT_SNIFFLES(SNIFFLES2.out.vcf.join(MINIMAP2.out.bam).join(MINIMAP2.out.bai), ch_reference, params.samplot_max_plots, 'sniffles2')

    // 7. Collect all files for MultiQC
    ch_multiqc_files = Channel.empty()
        .mix(MINIMAP2.out.stats.map { meta, stats -> stats }.collect().ifEmpty([]))

        // Add PEPPER stats if run
        .mix(ch_pepper_stats.stats.map { meta, stats -> stats }.collect().ifEmpty([]))
        .mix(ch_pepper_rtg.rtg_stats.map { meta, stats -> stats }.collect().ifEmpty([]))

        // Add CUTESV stats if run
        //.mix(ch_cutesv_stats.stats.map { meta, stats -> stats }.collect().ifEmpty([]))
        .mix(ch_cutesv_rtg.rtg_stats.map { meta, stats -> stats }.collect().ifEmpty([]))
        .mix(ch_cutesv_length_dist.stats.map { meta, stats -> stats }.collect().ifEmpty([]))
        .mix(ch_cutesv_length_dist.multiqc.ifEmpty([]))
        .mix(ch_cutesv_samplots.multiqc.ifEmpty([]))
        //.mix(ch_cutesv_stats.versions)

        // Add CLAIR3 stats if run
        .mix(ch_clair3_stats.stats.map { meta, stats -> stats }.collect().ifEmpty([]))
        .mix(ch_clair3_rtg.rtg_stats.map { meta, stats -> stats }.collect().ifEmpty([]))

        // Add SNIFFLES2 stats if run
        //.mix(ch_sniffles2_stats.stats.map { meta, stats -> stats }.collect().ifEmpty([]))
        .mix(ch_sniffles2_rtg.rtg_stats.map { meta, stats -> stats }.collect().ifEmpty([]))
        .mix(ch_sniffles2_length_dist.stats.map { meta, stats -> stats }.collect().ifEmpty([]))
        .mix(ch_sniffles2_length_dist.multiqc.ifEmpty([]))
        .mix(ch_sniffles2_samplots.multiqc.ifEmpty([]))
        //.mix(ch_sniffles2_stats.versions)

    // 8. Software versions
    // Each process writes one versions file per sample, all with the same
    // name, which makes MultiQC fail with a file name collision as soon as
    // there is more than one sample. De-duplicate by content and collate into
    // a single *_mqc_versions.yml, the file pattern MultiQC reads to build its
    // 'Software Versions' section.
    ch_versions = Channel.empty()
        .mix(MINIMAP2.out.versions)
        .mix(ch_pepper_stats.versions)
        .mix(ch_pepper_rtg.versions)
        .mix(PEPPER.out.versions)
        .mix(ch_cutesv_rtg.versions)
        .mix(CUTESV.out.versions)
        .mix(ch_clair3_stats.versions)
        .mix(ch_clair3_rtg.versions)
        .mix(CLAIR3.out.versions)
        .mix(ch_sniffles2_rtg.versions)
        .mix(SNIFFLES2.out.versions)

    ch_versions_collated = ch_versions
        .unique { it.text }
        .collectFile(name: 'software_versions_mqc_versions.yml', newLine: true, sort: true)

    ch_multiqc_files = ch_multiqc_files.mix(ch_versions_collated)

    // 9. Run MultiQC
    MULTIQC(
        ch_multiqc_files.collect(),
        Channel.value([]),  // multiqc_config (optional)
        Channel.value([]),  // extra_multiqc_config (optional)
        params.multiqc_logo ? Channel.fromPath(params.multiqc_logo) : Channel.value([]) // multiqc_logo (optional)
    )   
    
}

workflow.onComplete {
    println """
    Pipeline completed!
    Status:    ${workflow.success ? 'SUCCESS' : 'FAILED'}
    Work dir:  ${workflow.workDir}
    Results:   ${params.outdir}
    """
    
    // Only send when --email looks like an address (null, empty or the text "null" send nothing)
    if (params.email && params.email ==~ /\S+@\S+\.\S+/) {
        def status = workflow.success ? "SUCCESS" : "FAILED"
        def subject = "Long read pipeline status: ${status}"
        def body = """
        Pipeline execution summary
        ---------------------------
        Success           : ${workflow.success}
        Exit status       : ${workflow.exitStatus}
        Launch time       : ${workflow.start.format('dd-MMM-yyyy HH:mm:ss')}
        Ending time       : ${workflow.complete.format('dd-MMM-yyyy HH:mm:ss')} (duration: ${workflow.duration})
        Launch directory  : ${workflow.launchDir}
        Work directory    : ${workflow.workDir.toUriString()}
        Project directory : ${workflow.projectDir}
        Script ID         : ${workflow.scriptId ?: '-'}
        Workflow session  : ${workflow.sessionId}
        Nextflow run name : ${workflow.runName}
        Nextflow version  : ${workflow.nextflow.version}, build ${workflow.nextflow.build} (${workflow.nextflow.timestamp})
        
        Command:
        ${workflow.commandLine}

        ${workflow.success ? '' : """
        Errors:
        Error Message: ${workflow.errorMessage}
        Error Report : ${workflow.errorReport}
        """}
        """
        
        sendMail(to: params.email, subject: subject, body: body)
    }
}