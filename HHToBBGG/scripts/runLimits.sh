
text2workspace.py HHToBBGG_card_template.txt -m 125 

#limits
combine -M Asymptotic  -m 125 -d HHToBBGG_card_template.root

#significance
combine  -m 125 -M ProfileLikelihood --significance -d HHToBBGG_card_template.root -t -1 --expectSignal=1 --toysFreq




#Compute Limits
for s in {0..12}
do
   combine -M Asymptotic  -m 125 -d HHToBBGG_card_lum_$s.txt | tee HHToBBGG_limits_lum_$s.txt
done
for s in {0..12}
do
   combine -M Asymptotic  -m 125 -d HHToBBGG_card_pho_$s.txt | tee HHToBBGG_limits_pho_$s.txt
done
for s in {0..12}
do
   combine -M Asymptotic  -m 125 -d HHToBBGG_card_jet_$s.txt | tee HHToBBGG_limits_jet_$s.txt
done
for s in {0..5}
do
   combine -M Asymptotic  -m 125 -d HHToBBGG_card_endcapPhotonFakerate_$s.txt | tee HHToBBGG_limits_endcapPhotonFakerate_$s.txt
done
for s in {0..5}
do
   combine -M Asymptotic  -m 125 -d HHToBBGG_card_phoEff_$s.txt | tee HHToBBGG_limits_phoEff_$s.txt
done
for s in {0..5}
do
   combine -M Asymptotic  -m 125 -d HHToBBGG_card_btagEff_$s.txt | tee HHToBBGG_limits_btagEff_$s.txt
done



#Print Limit Results
for s in {0..12}
do
   cat HHToBBGG_limits_lum_$s.txt | grep "Median for expected limits" | awk '{print $5 ","}'
done
for s in {0..12}
do
   cat HHToBBGG_limits_pho_$s.txt | grep "Median for expected limits" | awk '{print $5 ","}'
done

for s in {0..12}
do
   cat HHToBBGG_limits_jet_$s.txt | grep "Median for expected limits" | awk '{print $5 ","}'
done

for s in {0..5}
do
   cat HHToBBGG_limits_endcapPhotonFakerate_$s.txt | grep "Median for expected limits" | awk '{print $5 ","}'
done
for s in {0..5}
do
   cat HHToBBGG_limits_phoEff_$s.txt | grep "Median for expected limits" | awk '{print $5 ","}'
done
for s in {0..5}
do
   cat HHToBBGG_limits_btagEff_$s.txt | grep "Median for expected limits" | awk '{print $5 ","}'
done

#Print Limit Results
for s in {0..12}
do
   cat HHToBBGG_limits_lum_$s.txt | grep "Sigma  for expected limits" | awk '{print $5 ","}'
done
for s in {0..12}
do
   cat HHToBBGG_limits_pho_$s.txt | grep "Sigma  for expected limits" | awk '{print $5 ","}'
done

for s in {0..12}
do
   cat HHToBBGG_limits_jet_$s.txt | grep "Sigma  for expected limits" | awk '{print $5 ","}'
done

for s in {0..5}
do
   cat HHToBBGG_limits_endcapPhotonFakerate_$s.txt | grep "Sigma  for expected limits" | awk '{print $5 ","}'
done
for s in {0..5}
do
   cat HHToBBGG_limits_phoEff_$s.txt | grep "Sigma  for expected limits" | awk '{print $5 ","}'
done
for s in {0..5}
do
   cat HHToBBGG_limits_btagEff_$s.txt | grep "Sigma  for expected limits" | awk '{print $5 ","}'
done
