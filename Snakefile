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
        "data-raw/HI_hanam.csv",
        "data-raw/Obj2_timecourse_200904_complete.xlsx",
        "data-raw/HI_long.csv",
        "data/data.R"
    output:
        "data/hi.csv",
        "data/hi-obj2.csv",
        "data/hi-rmh-hcw.csv",
        "data/hi-hanam.csv"
    shell:
        "Rscript data/data.R"

rule data_plot:
    input:
        ".deps-installed",
        "data/read_data.R",
        "data/hi.csv",
        "data/hi-obj2.csv",
        "data/hi-rmh-hcw.csv",
        "data/hi-hanam.csv",
        "data-plot/data-plot.R"
    output:
        directory("data-plot/indiv-hi"),
        directory("data-plot/indiv-hi-hanam"),
        directory("data-plot/indiv-hi-2"),
        directory("data-plot/indiv-hi-2-bwyears"),
        directory("data-plot/indiv-hi-rmh-hcw"),
        directory("data-plot/indiv-hi-alt"),
        directory("data-plot/indiv-hi-hanam-alt"),
        directory("data-plot/indiv-hi-2-alt"),
        directory("data-plot/indiv-hi-2-bwyears-alt"),
        directory("data-plot/indiv-hi-rmh-hcw-alt")
    shell:
        "Rscript data-plot/data-plot.R"

rule simple_diff:
    input:
        ".deps-installed",
        "data/hi.csv",
        "data/hi-rmh-hcw.csv",
        "data/read_data.R",
        "simple-diff/simple-diff.R"
    output:
        "simple-diff/simple-diff.png",
        "simple-diff/simple-diff-hanam.png",
        "simple-diff/simple-diff-hanam-d280.png",
        "simple-diff/simple-diff-rmh-hcw.png",
        "simple-diff/simple-diff.csv",
        "simple-diff/simple-diff-hanam.csv",
        "simple-diff/simple-diff-hanam-d280.csv",
        "simple-diff/simple-diff-rmh-hcw.csv"
    shell:
        "Rscript simple-diff/simple-diff.R"

rule drop_off:
    input:
        ".deps-installed",
        "data/hi-obj2.csv",
        "data/read_data.R",
        "drop-off/drop-off.R"
    output:
        "drop-off/spag.png",
        "drop-off/spag-even.png",
        "drop-off/spag-facets.png",
        "drop-off/vac-resp.png",
        "drop-off/vac-resp.csv"
    shell:
        "Rscript drop-off/drop-off.R"

rule all:
    input:
        rules.data.output,
        rules.data_plot.output,
        rules.simple_diff.output,
        rules.drop_off.output
