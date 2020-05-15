%THIS IS A DEMO FOR IMAGE ENHANCEMENT USING MEMETIC ALGORITHM EMPLOYING A 

load S.dat;
i = S;


%i=imresize(i,0.2);  %Demo run on resized image
%image = im2double(i);
image=i;
[row col]=size(image);
truelabels = image(:,1);
image = image(:,2:col);
[row col]=size(image);
tic();
MEMETICFS(image,truelabels);
toc();

     
      
    

