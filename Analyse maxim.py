# gemaakt door Eldo de Graaf, op 14/06 en 17/06
# dit is de code om data te analyseren, vervolgens is het de taak aan jullie (milan) om datapunten te verzamelen
# ik heb al 1 datapunt voorgedaan 

import matplotlib.pyplot as plt 

def analyse():
    
    
    #zet n op het aantal regels wat je hebt in je bestand
    n = 303
    # voer hier je gexporteerde file in van 03, zorg dat dit prc dezelfde coordinaten heeft als ha, of pas het aan zodat het overeen komt
    with open ('PiekO3.csv', 'r') as input_bestand_0III:
    
        
        
        l_intensiteit = []
        l_intensiteit_genormaliseerd = []
        l_pixelwaarde = []
        regelnummer = 0
        
        for regel in input_bestand_0III:
            regelnummer = regelnummer + 1
            
            if regelnummer >= 1 and regelnummer < (n+1):
                data_opgeknipt = regel.split(',')
                intensiteit_0 = (float(data_opgeknipt[1]))
                l_intensiteit.append(intensiteit_0)
        normalisatie = max(l_intensiteit)
        locatie_O3 = l_intensiteit.index(normalisatie)
        for j in range (0,n):
            l_intensiteit_genormaliseerd.append(l_intensiteit[j]/normalisatie)
            l_pixelwaarde.append(j)
        
    # voer hier je gexporteerde file in van ha, zorg dat dit prc dezelfde coordinaten heeft als 03, of pas het aan zodat het overeen komt
    with open ('PiekHa.csv', 'r') as input_bestand_Ha:
        
        l_intensiteit_Ha = []
        l_intensiteit_genormaliseerd_Ha = []
        regelnummer = 0
        
        for regel in input_bestand_Ha:
            regelnummer = regelnummer + 1
            
            if regelnummer >= 1 and regelnummer < (n+1):
                data_opgeknipt = regel.split(',')
                intensiteit_Ha = (float(data_opgeknipt[1]))
                l_intensiteit_Ha.append(intensiteit_Ha)
        normalisatie_Ha = max(l_intensiteit_Ha)
        lijst_i_Ha = []
        for i in l_intensiteit_Ha:
            lijst_i_Ha.append(i)
        lijst_i_Ha.sort()
        locatie_Ha = (l_intensiteit_Ha.index(lijst_i_Ha[-1]) + l_intensiteit_Ha.index(lijst_i_Ha[-2]) + l_intensiteit_Ha.index(lijst_i_Ha[-3]))/3
       
        for i in range (0,n):
            l_intensiteit_genormaliseerd_Ha.append(l_intensiteit_Ha[i]/normalisatie_Ha)
            
        
        # geeft verschil in pixels tussen beide pieken aan
        verschil_in_pixels = locatie_O3 - locatie_Ha
        schaal = 0.43
        verschil = verschil_in_pixels * schaal
        #print(verschil_in_pixels)
        print(f'verschil in gemiddelde afstand (2D) tussen pieken in arc seconds', verschil)
       
        plt.plot(l_pixelwaarde, l_intensiteit_genormaliseerd, 'g')
        plt.plot(l_pixelwaarde, l_intensiteit_genormaliseerd_Ha, 'r') 
        plt.xlabel('Pixelwaarde (mean verticale box)', fontsize = 10)
        plt.ylabel('intensiteit genormaliseerd', fontsize = 10)
        plt.legend(['0-III', 'Ha'], loc = 'upper right')
        plt.savefig('intensiteit.png')
        plt.show()
       
        
analyse()
