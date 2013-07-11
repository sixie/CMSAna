#############
# Signal
#############
cmsLs /store/user/sixie/BACON/V00-00-02/HHtoBBGG-8tev-START53_V7A/ | awk '{print $5}' >! CMSAna/HHToBBGG/catalog/HHtoBBGG-8tev-START53_V7A.txt
cmsLs /store/user/sixie/BACON/V00-00-02/HHtoBBGG-14tev-START53_V7A/ | awk '{print $5}' >! CMSAna/HHToBBGG/catalog/HHtoBBGG-14tev-START53_V7A.txt

#############
# Bkg
#############
cmsLs /store/user/sixie/BACON/V00-00-02/ttHgg-125-START53_V7A/ | awk '{print $5}' >! CMSAna/HHToBBGG/catalog/ttHgg-125-START53_V7A.txt
cmsLs /store/user/sixie/BACON/V00-00-02/ZHgg-125-START53_V7A/ | awk '{print $5}' >! CMSAna/HHToBBGG/catalog/ZHgg-125-START53_V7A.txt
cmsLs /store/user/sixie/BACON/V00-00-02/ggHgg-125-START53_V7A/ | awk '{print $5}' >! CMSAna/HHToBBGG/catalog/ggHgg-125-START53_V7A.txt

