rule install_deps:
    input:
        "renv.lock"
    output:
        touch(".deps-installed")
    shell:
        """Rscript -e 'renv::restore()'"""

rule data:
    input:
        ".deps-installed",
        "data-raw/HI.xlsx",
        "data-raw/Obj1_Sample_info.xlsx",
        "data-raw/VirusDiln_BT.xlsx",
        "data-raw/Viruses.xlsx",
        "data-raw/HI-annette-extra.csv",
        "data-raw/Obj2_timecource.xlsx",
        "data-raw/HI_long.csv",
        "data/data.R"
    output:
        "data/hi.csv",
        "data/hi-rmh-hcw.csv",
        "data/hi-annette-extra.csv"
    shell:
        "Rscript data/data.R"

rule data_plot:
    input:
        ".deps-installed",
        "data/hi.csv",
        "data/read_data.R",
        "data-plot/data-plot.R"
    output:
        directory("data-plot/indiv-hi"),
        directory("data-plot/indiv-hi-annette-extra"),
        directory("data-plot/indiv-hi-alt"),
        directory("data-plot/indiv-hi-alt-annette-extra")
    shell:
        "Rscript data-plot/data-plot.R"

rule simple_diff:
    input:
        ".deps-installed",
        "data/hi.csv",
        "data/read_data.R",
        "simple-diff/simple-diff.R"
    output:
        "simple-diff/simple-diff.pdf"
    shell:
        "Rscript simple-diff/simple-diff.R"

rule all:
    input:
        rules.data.output,
        rules.data_plot.output,
        rules.simple_diff.output
