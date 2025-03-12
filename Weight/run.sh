cd /Users/xiongzheng/software/B4/B4c/Weight
# cp /Users/xiongzheng/software/B4/B4c/Root/Electron_PowerLaw.root   /Users/xiongzheng/software/B4/B4c/Weight/Electron_PowerLaw.root  
# cp /Users/xiongzheng/software/B4/B4c/Root/Proton_PowerLaw.root     /Users/xiongzheng/software/B4/B4c/Weight/Proton_PowerLaw.root    
# cp /Users/xiongzheng/software/B4/B4c/Root/Deuteron_PowerLaw.root   /Users/xiongzheng/software/B4/B4c/Weight/Deuteron_PowerLaw.root  
# cp /Users/xiongzheng/software/B4/B4c/Root/Helium3_PowerLaw.root    /Users/xiongzheng/software/B4/B4c/Weight/Helium3_PowerLaw.root   
# cp /Users/xiongzheng/software/B4/B4c/Root/Helium4_PowerLaw.root    /Users/xiongzheng/software/B4/B4c/Weight/Helium4_PowerLaw.root   
cp /Users/xiongzheng/software/B4/B4c/Root/Carbon_PowerLaw.root     /Users/xiongzheng/software/B4/B4c/Weight/Carbon_PowerLaw.root    


exe=/Users/xiongzheng/software/B4/B4c/Weight/Normalized
filepath=/Users/xiongzheng/software/B4/B4c/Weight/
Tag=(Carbon) # Carbon Helium3 Proton Deuteron Electron Helium4
for i in "${!Tag[@]}"; 
do
    echo "${exe} ${Tag[$i]}_PowerLaw.root"
    time ${exe}  ${Tag[$i]}_PowerLaw.root 
done

