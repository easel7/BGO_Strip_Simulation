cd /Users/xiongzheng/software/B4/B4e/Weight
# cp /Users/xiongzheng/software/B4/B4e/Root/Electron_PowerLaw.root   /Users/xiongzheng/software/B4/B4e/Weight/Electron_PowerLaw.root  
# cp /Users/xiongzheng/software/B4/B4e/Root/Proton_PowerLaw.root     /Users/xiongzheng/software/B4/B4e/Weight/Proton_PowerLaw.root    
# cp /Users/xiongzheng/software/B4/B4e/Root/Deuteron_PowerLaw.root   /Users/xiongzheng/software/B4/B4e/Weight/Deuteron_PowerLaw.root  
# cp /Users/xiongzheng/software/B4/B4e/Root/Helium3_PowerLaw.root    /Users/xiongzheng/software/B4/B4e/Weight/Helium3_PowerLaw.root   
# cp /Users/xiongzheng/software/B4/B4e/Root/Helium4_PowerLaw.root    /Users/xiongzheng/software/B4/B4e/Weight/Helium4_PowerLaw.root   
cp /Users/xiongzheng/software/B4/B4e/Root/Carbon_PowerLaw.root     /Users/xiongzheng/software/B4/B4e/Weight/Carbon_PowerLaw.root    



exe=/Users/xiongzheng/software/B4/B4e/Weight/Normalized
filepath=/Users/xiongzheng/software/B4/B4e/Weight/
Tag=(Carbon ) # Carbon Helium3 Proton Deuteron Electron Helium4 Helium3 Proton Deuteron Electron Helium4
for i in "${!Tag[@]}"; 
do
    echo "${exe} ${Tag[$i]}_PowerLaw.root"
    time ${exe}  ${Tag[$i]}_PowerLaw.root 
done

