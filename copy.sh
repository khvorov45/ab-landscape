PROJ="/home/khvorova/Projects/ab-landscape"
SHARED="/home/khvorova/vidrlwhoflu/Research"

OBJ1="$SHARED/BAA 75D301-10-R-67835/Analyses/Objective 1"
OBJ2="$SHARED/BAA 75D301-10-R-67835/Analyses/Objective 2"

rm -rf "$OBJ2/Israel/indiv-hi-2" "$OBJ2/Israel/indiv-hi-2-alt"
cp -R $PROJ/data-plot/indiv-hi-2 $PROJ/data-plot/indiv-hi-2-alt "$OBJ2/Israel"

rm -rf "$OBJ2/Israel/indiv-hi-2-bwyears" "$OBJ2/Israel/indiv-hi-2-bwyears-alt"
cp -R $PROJ/data-plot/indiv-hi-2-bwyears $PROJ/data-plot/indiv-hi-2-bwyears-alt "$OBJ2/Israel"

rm -rf "$OBJ1/Israel/indiv-hi" "$OBJ1/Israel/indiv-hi-alt"
cp -R $PROJ/data-plot/indiv-hi $PROJ/data-plot/indiv-hi-alt "$OBJ1/Israel"

rm -rf "$SHARED/Ha Nam/landscape-plots/indiv-hi-hanam" "$SHARED/Ha Nam/landscape-plots/indiv-hi-hanam-alt"
cp -R $PROJ/data-plot/indiv-hi-hanam $PROJ/data-plot/indiv-hi-hanam-alt "$SHARED/Ha Nam/landscape-plots"

rm -rf "$SHARED/2016.003 RMH HCW Study/Analyses/landscape-plots/indiv-hi-rmh-hcw" "$SHARED/2016.003 RMH HCW Study/Analyses/landscape-plots/indiv-hi-rmh-hcw-alt"
cp -R $PROJ/data-plot/indiv-hi-rmh-hcw $PROJ/data-plot/indiv-hi-rmh-hcw-alt "$SHARED/2016.003 RMH HCW Study/Analyses/landscape-plots"


cp $PROJ/simple-diff/simple-diff-rmh-hcw.png "$SHARED/2016.003 RMH HCW Study/Analyses/landscape-plots/simple-diff-rmh-hcw.png"

cp $PROJ/simple-diff/simple-diff-hanam.png "$SHARED/Ha Nam/landscape-plots/simple-diff-hanam.png"
cp $PROJ/simple-diff/simple-diff-hanam-d280.png "$SHARED/Ha Nam/landscape-plots/simple-diff-hanam-d280.png"

cp $PROJ/simple-diff/simple-diff.png "$OBJ1/Israel/simple-diff.png"


rm -rf "$OBJ2/drop-off/*"
cp "$PROJ/drop-off/spag.png" "$PROJ/drop-off/spag-even.png" "$PROJ/drop-off/spag-facets.png" "$PROJ/drop-off/vac-resp.png" "$PROJ/drop-off/vac-resp.csv" "$PROJ/drop-off/vac-resp-scatter.png" "$PROJ/drop-off/vac-resp-cat.png" "$PROJ/drop-off/vac-resp-vacc.png" "$PROJ/drop-off/vac-resp-age.png" "$OBJ2/drop-off"

rm -rf /home/khvorova/vidrlwhoflu/Research/ab-landscape
cp -R /home/khvorova/Projects/ab-landscape /home/khvorova/vidrlwhoflu/Research/ab-landscape
