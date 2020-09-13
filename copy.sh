PROJ="/home/khvorova/Projects/ab-landscape"
SHARED="/home/khvorova/vidrlwhoflu/Research"

OBJ1="$SHARED/BAA 75D301-10-R-67835/Analyses/Objective 1"
OBJ2="$SHARED/BAA 75D301-10-R-67835/Analyses/Objective 2"

rm -rf "/home/khvorova/vidrlwhoflu/Research/BAA 75D301-10-R-67835/Analyses/Objective 2/Israel/indiv-hi-2" "/home/khvorova/vidrlwhoflu/Research/BAA 75D301-10-R-67835/Analyses/Objective 2/Israel/indiv-hi-2-alt"
cp -R /home/khvorova/Projects/ab-landscape/data-plot/indiv-hi-2 /home/khvorova/Projects/ab-landscape/data-plot/indiv-hi-2-alt "/home/khvorova/vidrlwhoflu/Research/BAA 75D301-10-R-67835/Analyses/Objective 2/Israel"

rm -rf "/home/khvorova/vidrlwhoflu/Research/BAA 75D301-10-R-67835/Analyses/Objective 2/Israel/indiv-hi-2-bwyears" "/home/khvorova/vidrlwhoflu/Research/BAA 75D301-10-R-67835/Analyses/Objective 2/Israel/indiv-hi-2-bwyears-alt"
cp -R /home/khvorova/Projects/ab-landscape/data-plot/indiv-hi-2-bwyears /home/khvorova/Projects/ab-landscape/data-plot/indiv-hi-2-bwyears-alt "/home/khvorova/vidrlwhoflu/Research/BAA 75D301-10-R-67835/Analyses/Objective 2/Israel"

rm -rf "/home/khvorova/vidrlwhoflu/Research/BAA 75D301-10-R-67835/Analyses/Objective 1/Israel/indiv-hi" "/home/khvorova/vidrlwhoflu/Research/BAA 75D301-10-R-67835/Analyses/Objective 1/Israel/indiv-hi-alt"
cp -R /home/khvorova/Projects/ab-landscape/data-plot/indiv-hi /home/khvorova/Projects/ab-landscape/data-plot/indiv-hi-alt "/home/khvorova/vidrlwhoflu/Research/BAA 75D301-10-R-67835/Analyses/Objective 1/Israel"

rm -rf "/home/khvorova/vidrlwhoflu/Research/Ha Nam/landscape-plots/indiv-hi-hanam" "/home/khvorova/vidrlwhoflu/Research/Ha Nam/landscape-plots/indiv-hi-hanam-alt"
cp -R /home/khvorova/Projects/ab-landscape/data-plot/indiv-hi-hanam /home/khvorova/Projects/ab-landscape/data-plot/indiv-hi-hanam-alt "/home/khvorova/vidrlwhoflu/Research/Ha Nam/landscape-plots"

rm -rf "/home/khvorova/vidrlwhoflu/Research/2016.003 RMH HCW Study/Analyses/landscape-plots/indiv-hi-rmh-hcw" "/home/khvorova/vidrlwhoflu/Research/2016.003 RMH HCW Study/Analyses/landscape-plots/indiv-hi-rmh-hcw-alt"
cp -R /home/khvorova/Projects/ab-landscape/data-plot/indiv-hi-rmh-hcw /home/khvorova/Projects/ab-landscape/data-plot/indiv-hi-rmh-hcw-alt "/home/khvorova/vidrlwhoflu/Research/2016.003 RMH HCW Study/Analyses/landscape-plots"


cp /home/khvorova/Projects/ab-landscape/simple-diff/simple-diff-rmh-hcw.png "/home/khvorova/vidrlwhoflu/Research/2016.003 RMH HCW Study/Analyses/landscape-plots/simple-diff-rmh-hcw.png"

cp /home/khvorova/Projects/ab-landscape/simple-diff/simple-diff-hanam.png "/home/khvorova/vidrlwhoflu/Research/Ha Nam/landscape-plots/simple-diff-hanam.png"
cp /home/khvorova/Projects/ab-landscape/simple-diff/simple-diff-hanam-d280.png "/home/khvorova/vidrlwhoflu/Research/Ha Nam/landscape-plots/simple-diff-hanam-d280.png"

cp /home/khvorova/Projects/ab-landscape/simple-diff/simple-diff.png "/home/khvorova/vidrlwhoflu/Research/BAA 75D301-10-R-67835/Analyses/Objective 1/Israel/simple-diff.png"


rm -rf "$OBJ2/drop-off/*"
cp "$PROJ/drop-off/spag.png" "$PROJ/drop-off/spag-even.png" "$PROJ/drop-off/spag-facets.png" "$PROJ/drop-off/vac-resp.png" "$PROJ/drop-off/vac-resp.csv" "$PROJ/drop-off/vac-resp-scatter.png" "$PROJ/drop-off/vac-resp-cat.png" "$PROJ/drop-off/vac-resp-vacc.png" "$PROJ/drop-off/vac-resp-age.png" "$OBJ2/drop-off"

rm -rf /home/khvorova/vidrlwhoflu/Research/ab-landscape
cp -R /home/khvorova/Projects/ab-landscape /home/khvorova/vidrlwhoflu/Research/ab-landscape
