pushd /fs/ncanda-share/cases

for CASE in NCANDA_S0*/standard/*; do 
  FILE=$CASE/restingstate/rs_specific_processing/qa/art/art.bold_noIntenCorr_4d_outliers.txt; 
  if [ -e $FILE ]; then 
    echo "`echo $CASE | tr '/' ','`,`wc -l $FILE| cut -d' ' -f1`"; 
  fi; 
done 

