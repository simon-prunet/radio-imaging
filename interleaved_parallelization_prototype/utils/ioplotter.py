import sys
import matplotlib.pyplot as plt
import numpy as np

# Vérifier si le fichier d'entrée a été fourni                                                                                                                                                                                         
if len(sys.argv) < 3:
    print("Usage: python plot_read_write.py input_file output_file")
    sys.exit(1)

# Charger les données depuis le fichier fourni en argument, en sélectionnant les colonnes pertinentes (temps, kB_rd/s, kB_wr/s)                                                                                                        
input_file = sys.argv[1]
output_file = sys.argv[2]

data = np.loadtxt(input_file, usecols=(0, 3, 4))

# Extraire les colonnes                                                                                                                                                                                                                
time = data[:, 0] - data[0, 0]  # Temps                                                                                                                                                                                                             
kB_rd_s = data[:, 1]  # kB read/s                                                                                                                                                                                                      
kB_wr_s = data[:, 2]  # kB write/s                                                                                                                                                                                                     
# Largeur des barres                                                                                                                                                                                                                   
bar_width = 0.4
print(np.sum(kB_rd_s))
print(np.sum(kB_wr_s))
# Créer un histogramme pour kB read/s et kB write/s                                                                                                                                                                                    
plt.figure(figsize=(12, 6))

# Histogramme pour kB read/s                                                                                                                                                                                                           
plt.bar(time - bar_width / 2, kB_rd_s / 1024, width=bar_width, label='MB_read/s', color='darkblue', edgecolor='darkblue', alpha=0.8)

# Histogramme pour kB write/s                                                                                                                                                                                                          
plt.bar(time + bar_width / 2, kB_wr_s / 1024, width=bar_width, label='MB_wr/s', color='orange', edgecolor='orange', alpha=0.8)


# Configurer les axes                                                                                                                                                                                                                  
plt.xlabel('Time (s)')
plt.ylabel('MB/s')
plt.title('I/O Rate over Time (MB/s)')
plt.legend()
plt.yscale('log')

# Ajouter une grille                                                                                                                                                                                                                   
plt.grid(axis='y', linestyle='--', alpha=0.7)

# Imposer l'échelle des ordonnées de 0 à 3000 MB/s                                                                                                                                                                                     
#plt.ylim(0, 100)

# Fixer la limite des abscisses à 1200 secondes                                                                                                                                                                                        
#plt.xlim(0, 1200)                                                                                                                                                                                                                     
#plt.xlim(left=-50, right=1200)
# Configurer les graduations de l'axe des ordonnées                                                                                                                                                                                    
#plt.yticks(np.arange(0, 3000 + 1, 500),  # Graduations de 500 en 500                                                                                                                                                                   
#           labels=np.arange(0, 3000 + 1, 500).astype(int))

# Ajuster la mise en page                                                                                                                                                                                                              
plt.tight_layout()
plt.savefig(output_file)

# Afficher le graphique                                                                                                                                                                                                                
plt.show()