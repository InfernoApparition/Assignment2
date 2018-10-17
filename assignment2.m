%Whole code has been written from scratch. No references were taken while
%writing the code.

global kernel kernelmatrix img padmatrix K
[address,No_img_selected] = imgetfile(); %function return two variable "image address", "image selected or not" 
kernel = imread(address);
if No_img_selected
    msgbox(springf('Please select a kernel'));
    return
end
kernel = im2double(imread(address));

kernel1 = kernel(1:81 , 1:70);
kernel1sum=sum(kernel1(:));
kernelmatrix=kernel1./kernel1sum;
 
[address,No_img_selected] = imgetfile(); %function return two variable "image address", "image selected or not" 
img = im2double(imread(address));
if No_img_selected
    msgbox(springf('Please select an image'));
    return
end
imshow(img);
title('Image');

rchannel=img(:,:,1);
gchannel=img(:,:,2);
bchannel=img(:,:,3);
[rows,columns]=size(rchannel);


%Full Inverse Filter
padmatrix = zeros(rows,columns);
for i=1:81
    for j=1:70
        padmatrix(i+floor((rows-80)/2),j+floor((columns-71)/2))=kernelmatrix(i,j);
    end
end

fftpadmatrix=fft2(padmatrix);

rblur = ifft2(fft2(rchannel).*fftpadmatrix);
gblur = ifft2(fft2(gchannel).*fftpadmatrix);
bblur = ifft2(fft2(bchannel).*fftpadmatrix);
rbluro=rblur;
gbluro=gblur;
bbluro=bblur;
for i=1:rows
    for j=1:columns
        if(i<floor(rows/2) && j>floor(columns/2)) rblur(i+floor(rows/2),j-floor(columns/2))=rbluro(i,j); gblur(i+floor(rows/2),j-floor(columns/2))=gbluro(i,j); bblur(i+floor(rows/2),j-floor(columns/2))=bbluro(i,j);end
        if(i<floor(rows/2) && j<floor(columns/2)) rblur(i+floor(rows/2),j+floor(columns/2))=rbluro(i,j); gblur(i+floor(rows/2),j+floor(columns/2))=gbluro(i,j); bblur(i+floor(rows/2),j+floor(columns/2))=bbluro(i,j);end
        if(i>floor(rows/2) && j>floor(columns/2)) rblur(i-floor(rows/2),j-floor(columns/2))=rbluro(i,j); gblur(i-floor(rows/2),j-floor(columns/2))=gbluro(i,j); bblur(i-floor(rows/2),j-floor(columns/2))=bbluro(i,j);end
        if(i>floor(rows/2) && j<floor(columns/2)) rblur(i-floor(rows/2),j+floor(columns/2))=rbluro(i,j); gblur(i-floor(rows/2),j+floor(columns/2))=gbluro(i,j); bblur(i-floor(rows/2),j+floor(columns/2))=bbluro(i,j);end
    end
end
blurimg(:,:,1)=rblur;
blurimg(:,:,2)=gblur;
blurimg(:,:,3)=bblur;
blurimg1(:,:,1)=rbluro;
blurimg1(:,:,2)=gbluro;
blurimg1(:,:,3)=bbluro;
figure, imshow(blurimg);
title('Blurred Image');

ftrblur1=fft2(rblur);
ftgblur1=fft2(gblur);
ftbblur1=fft2(bblur);

%FT of blurred image
ftrblur=fft2(rbluro);
ftgblur=fft2(gbluro);
ftbblur=fft2(bbluro);

for i=1:rows
    for j=1:columns
        if(fftpadmatrix(i,j)==0)
            fftpadmatrix(i,j)=1;
        end
    end
end

%divide by kernel

finalfftr=ifft2(ftrblur./fftpadmatrix);
finalfftg=ifft2(ftgblur./fftpadmatrix);
finalfftb=ifft2(ftbblur./fftpadmatrix);

finalimg(:,:,1)=finalfftr;
finalimg(:,:,2)=finalfftg;
finalimg(:,:,3)=finalfftb;

figure, imshow(finalimg);
title('Full Inverse Image');

%Truncated Inverse Filter
trunr = ftrblur./fftpadmatrix;
trung = ftgblur./fftpadmatrix;
trunb = ftbblur./fftpadmatrix;

prompt = 'r Value: ';%Input r Value
input = inputdlg(prompt);
r = str2num(input{1});

for i=1:rows
    for j=1:columns
        if(i<floor(rows/2) && j>floor(columns/2))
            if (i^2+(columns-j)^2>r^2) trunr(i,j)=1; trung(i,j)=1; trunb(i,j)=1; 
            end
        end
        if(i<floor(rows/2) && j<floor(columns/2))
            if (i^2+j^2>r^2) trunr(i,j)=1; trung(i,j)=1; trunb(i,j)=1;
            end
        end
        if(i>floor(rows/2) && j>floor(columns/2))
            if ((rows-i)^2+(columns-j)^2>r^2) trunr(i,j)=1; trung(i,j)=1; trunb(i,j)=1;
            end
        end
        if(i>floor(rows/2) && j<floor(columns/2))
            if ((rows-i)^2+j^2>r^2) trunr(i,j)=1; trung(i,j)=1; trunb(i,j)=1;
            end
        end
    end
end

truniftr=ifft2(trunr);
truniftg=ifft2(trung);
truniftb=ifft2(trunb);

trunf(:,:,1)=truniftr;
trunf(:,:,2)=truniftg;
trunf(:,:,3)=truniftb;

figure, imshow(trunf);
title('Truncated Inverse Image');
[peaksnrTIF, snr] = psnr(trunf, img);
ReTIFimg = real(trunf);
TIFssimval = ssim(ReTIFimg, img);
fprintf('\n The Peak-SNR value after trauncated inverse filter is %0.4f', peaksnrTIF);
fprintf('\n The SSIM after truncated filter is %0.4f', TIFssimval);

%Weiner
Hsquared=abs(fftpadmatrix).^2;
prompt = 'K Value: ';%Input K Value
input = inputdlg(prompt);
k = str2num(input{1});
K(1:rows,1:columns) = k;
comjH=conj(fftpadmatrix);

%{
fimgr= (Hsquared.*rbluro)./(fftpadmatrix.*(K+Hsquared));
fimgg= (Hsquared.*gbluro)./(fftpadmatrix.*(K+Hsquared));
fimgb= (Hsquared.*bbluro)./(fftpadmatrix.*(K+Hsquared));
%}

fimgr= ifft2((comjH.*ftrblur)./(K+Hsquared));
fimgg= ifft2((comjH.*ftgblur)./(K+Hsquared));
fimgb= ifft2((comjH.*ftbblur)./(K+Hsquared));

fimgf(:,:,1)=fimgr;
fimgf(:,:,2)=fimgg;
fimgf(:,:,3)=fimgb;

figure, imshow(fimgf);
title('Weiner Image');
[peaksnrWF, snr] = psnr(fimgf, img);
ReWFimg = real(fimgf);
WFssimval = ssim(ReWFimg, img);
fprintf('\n The Peak-SNR value after Weiner filter is %0.4f', peaksnrWF);
fprintf('\n The SSIM after Weiner filter is %0.4f', WFssimval);

%Constrained Least Square Filter
prompt = 'Give gamma Value: ';%Input Gamma Value
input = inputdlg(prompt);
gamma = str2num(input{1});
%for gamma = 0:0.1:1
pmatrix=[0 -1 0 ; -1 4 -1 ; 0 -1 0];
ppadmatrix=zeros(rows,columns);
for i=1:3
    for j=1:3
        ppadmatrix(i+floor((rows-3)/2),j+floor((columns-3)/2))=pmatrix(i,j);
    end
end

fftppadmatrix=(abs(fft2(ppadmatrix)).^2).*gamma;

f3imgr = ifft2((comjH.*ftrblur)./(Hsquared+fftppadmatrix));
f3imgg = ifft2((comjH.*ftgblur)./(Hsquared+fftppadmatrix));
f3imgb = ifft2((comjH.*ftbblur)./(Hsquared+fftppadmatrix));

f3imgf(:,:,1)= f3imgr;
f3imgf(:,:,2)= f3imgg;
f3imgf(:,:,3)= f3imgb;

figure, imshow(f3imgf);
title('CLSR Image');
[peaksnrCLSR, snr] = psnr(f3imgf, img);
%PSNR_CLSR(i) = peaksnr_CSLR;
ReCLSRimg = real(f3imgf); 
CLSRssimval = ssim(ReCLSRimg, img);
fprintf('\n The Peak-SNR value for CLSR filter is %0.4f', peaksnrCLSR);
fprintf('\n The SSIM for CSLR filter is %0.4f', CLSRssimval);