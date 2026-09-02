function [line]=get_the_line_profile(star,starname,type,s,files)

%Build filename
file_name=[files.list{s} '.fit'];

%Load data
A=fitsread(file_name); %load data

%check_the_file(A,s) %check quality of file

%Check image size
row_range=1960:2039;
if size(A,1)<max(row_range)
    error(['There is a size issue with file: ' file_name])
end

%creat matrix with 80 rows to average around 2000 (middle)
y=A(row_range,:);

%AVERAGE ALONG COLOUMNS
avrg=mean(y,1,'omitnan');%average through the coloums

%invert array and make for y
y=avrg';

%index the double array for use later
y=double(y(:,1));

%get data to lie on the floor
y=y-min(y);

%normalise the data
maximum_line_value=max(y);
if ~isfinite(maximum_line_value) || maximum_line_value<=0
    line=zeros(size(y));
else
    line=y/maximum_line_value;
end

end