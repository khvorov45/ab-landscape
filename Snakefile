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
        "data/data.R"
    output:
        "data/hi.csv"
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
        directory("data-plot/indiv-hi-alt")
    shell:
        "Rscript data-plot/data-plot.R"

rule all:
    input:
        rules.data.output,
        rules.data_plot.output
