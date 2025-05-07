exe=/Users/xiongzheng/software/B4/B4e/build/exampleB4e
Tag=(Proton Deuteron) #Proton Deuteron Electron Helium4 Carbon Helium3  Deuteron

for i in "${!Tag[@]}"; 
do
    echo  /Users/xiongzheng/software/B4/B4e/Macro/${Tag[$i]}_PowerLaw.mac
    ${exe} -m  /Users/xiongzheng/software/B4/B4e/Macro/${Tag[$i]}_PowerLaw.mac
done