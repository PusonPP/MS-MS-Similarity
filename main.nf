nextflow.enable.dsl = 2

params.raw = "${projectDir}/data/raw"
params.manifest = "${projectDir}/data/manifest.json"
params.settings = "${projectDir}/params.json"
params.outdir = "${projectDir}/results"

process PREPROCESS_SPECTRA {
    publishDir params.outdir, mode: 'copy'
    input:
    path raw
    path manifest
    path settings
    path code
    output:
    path 'preprocess', emit: data
    script:
    """
    python ${code}/preprocess.py ${raw} ${manifest} ${settings} preprocess
    """
}

process SCORE_SPECTRA {
    publishDir "${params.outdir}/scores", mode: 'copy'
    input:
    path clean
    path settings
    path code
    output:
    path 'spectral_scores.csv'
    script:
    """
    python ${code}/score.py ${clean}/spectra.json ${settings} spectral_scores.csv
    """
}

process COMPARE_STRUCTURES {
    publishDir params.outdir, mode: 'copy'
    input:
    path clean
    path scores
    path settings
    path code
    output:
    path 'comparison'
    script:
    """
    python ${code}/compare.py ${clean}/spectra.json ${scores} ${settings} comparison
    """
}

process BUILD_REPORT {
    publishDir params.outdir, mode: 'copy'
    input:
    path clean
    path comparison
    path code
    output:
    path 'report'
    script:
    """
    python ${code}/report.py ${clean} ${comparison} report
    """
}

workflow {
    raw = file(params.raw, checkIfExists: true)
    manifest = file(params.manifest, checkIfExists: true)
    settings = file(params.settings, checkIfExists: true)
    code = file("${projectDir}/scripts", checkIfExists: true)
    PREPROCESS_SPECTRA(raw, manifest, settings, code)
    SCORE_SPECTRA(PREPROCESS_SPECTRA.out.data, settings, code)
    COMPARE_STRUCTURES(PREPROCESS_SPECTRA.out.data, SCORE_SPECTRA.out, settings, code)
    BUILD_REPORT(PREPROCESS_SPECTRA.out.data, COMPARE_STRUCTURES.out, code)
}
