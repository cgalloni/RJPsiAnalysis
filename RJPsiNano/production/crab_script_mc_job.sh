echo "================= CMSRUN starting ===================="
#date=$(date '+%Y%b%d')
cmsRun -j FrameworkJobReport.xml -p PSet.py 
#cmsRun -j FrameworkJobReport.xml -p ../test/run_nano_jpsi_cfg.py
echo "================= CMSRUN finished ===================="
#right name for the input file
echo "==================puReweight starting ================"
#echo ${date} 

python  puReweight_2016.py --tag 2023Oct09
# output file slightly different name
#rm input file
echo "================= puReweight finished ================"


