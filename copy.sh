PROJ="/home/khvorova/Projects/ab-landscape"
SHARED="/home/khvorova/vidrlwhoflu/Research"

OBJ1="$SHARED/BAA 75D301-10-R-67835/Analyses/Objective 1"
OBJ2="$SHARED/BAA 75D301-10-R-67835/Analyses/Objective 2"

RMH="$SHARED/2016.003 RMH HCW Study"

HANAM="$SHARED/Ha Nam"

# Objective 1 -----------------------------------------------------------------

# HI plots
rm -rf "$OBJ1/Israel/indiv-hi" \
    "$OBJ1/Israel/indiv-hi-alt" \
    "$OBJ1/Israel/indiv-hi-contour"
cp -R "$PROJ/data-plot/indiv-hi" \
    "$PROJ/data-plot/indiv-hi-alt" \
    "$PROJ/data-plot/indiv-hi-contour" \
    "$OBJ1/Israel"

# Differences
cp $PROJ/simple-diff/simple-diff.png "$OBJ1/Israel/simple-diff.png"

# Objective 2 -----------------------------------------------------------------

# HI plots
rm -rf "$OBJ2/Israel/indiv-hi-2" \
    "$OBJ2/Israel/indiv-hi-2-alt" \
    "$OBJ2/Israel/indiv-hi-2-bwyears" \
    "$OBJ2/Israel/indiv-hi-2-bwyears-alt"
cp -R "$PROJ/data-plot/indiv-hi-2" \
    "$PROJ/data-plot/indiv-hi-2-alt" \
    "$PROJ/data-plot/indiv-hi-2-bwyears" \
    "$PROJ/data-plot/indiv-hi-2-bwyears-alt" \
    "$OBJ2/Israel"

# Drop-off evaluation
rm -rf "$OBJ2/drop-off/*"
cp "$PROJ/drop-off/spag.png" \
    "$PROJ/drop-off/spag-even.png" \
    "$PROJ/drop-off/spag-facets.png" \
    "$PROJ/drop-off/vac-resp.png" \
    "$PROJ/drop-off/vac-resp.csv" \
    "$PROJ/drop-off/vac-resp-scatter.png" \
    "$PROJ/drop-off/vac-resp-cat.png" \
    "$PROJ/drop-off/vac-resp-vacc.png" \
    "$PROJ/drop-off/vac-resp-age.png" \
    "$OBJ2/drop-off"

# Hanam -----------------------------------------------------------------------

# HI plots
rm -rf "$HANAM/landscape-plots/indiv-hi-hanam" \
    "$HANAM/landscape-plots/indiv-hi-hanam-alt"
cp -R "$PROJ/data-plot/indiv-hi-hanam" \
    "$PROJ/data-plot/indiv-hi-hanam-alt" \
    "$HANAM/landscape-plots"

# Differences
cp "$PROJ/simple-diff/simple-diff-hanam.png" \
    "$HANAM/landscape-plots/simple-diff-hanam.png"
cp "$PROJ/simple-diff/simple-diff-hanam-d280.png" \
    "$HANAM/landscape-plots/simple-diff-hanam-d280.png"

# RMH -------------------------------------------------------------------------

# HI plots
rm -rf "$RMH/Analyses/landscape-plots/indiv-hi-rmh-hcw" \
    "$RMH/Analyses/landscape-plots/indiv-hi-rmh-hcw-alt"
cp -R "$PROJ/data-plot/indiv-hi-rmh-hcw" \
    "$PROJ/data-plot/indiv-hi-rmh-hcw-alt" \
    "$RMH/Analyses/landscape-plots"

# Differences
cp "$PROJ/simple-diff/simple-diff-rmh-hcw.png" \
    "$RMH/Analyses/landscape-plots/simple-diff-rmh-hcw.png"

# The whole project -----------------------------------------------------------

rm -rf "$SHARED/ab-landscape"
cp -R $PROJ "$SHARED/ab-landscape"
