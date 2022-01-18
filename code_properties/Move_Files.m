raw_dir = '~/Data/Matteo_Early/'; 
glm_dir = '~/Data/Matteo_GLMEarly/'; 
process_dir = '/home/fstella/Dropbox/Projects_NIJ/Matteo/Nijmegen_LinearEarly/';

mkdir([process_dir 'Cluster_Info/'])

for animal = [2 3 5]

copyfile([raw_dir 'Animal' num2str(animal) '/Cluster*'] , [glm_dir 'Animal' num2str(animal)])
copyfile([raw_dir 'Animal' num2str(animal) '/Cluster*'] , [glm_dir 'Cluster_Info/'])


end