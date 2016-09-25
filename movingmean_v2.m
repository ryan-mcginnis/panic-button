function result=movingmean_v2(data,window,dim)
%Calculates the centered moving average of an n-dimensional matrix in any direction. 
%   result=movingmean(data,window,dim,option)

%   Inputs: 
%   1)  data = The matrix to be averaged. 
%   2)  window = The window size.  This works best as an odd number.  If an even 
%       number is entered then it is rounded down to the next odd number.  
%   3)  dim = The dimension in which you would like do the moving average.
%       This is an optional input.  To leave blank use [] place holder in
%       function call.  Defaults to 1.
%   4)  option = which solution algorithm to use.  The default option works
%       best in most situations, but option 2 works better for wide
%       matrices (i.e. 1000000 x 1 or 10 x 1000 x 1000) when solved in the
%       shorter dimension.  Data size where option 2 is more efficient will 
%       vary from computer to computre.This is an optional input.  To leave  
%       blank use [] place holder in function call.  Defaults to 1.
% 
%   Example:  
%   Calculate column moving average of 10000 x 10 matrix with a window size 
%   of 5 in the 1st dimension using algorithm option 1.
%   d=rand(10000,10);
%   dd=movingmean(d,5,1,1);
%           or
%   dd=movingmean(d,5,[],1);
%           or
%   dd=movingmean(d,5,1,[]);
%           or
%   dd=movingmean(d,5,1);
%           or
%   dd=movingmean(d,5,[],[]);
%           or
%   dd=movingmean(d,5);
%
%   Moving mean for each element uses data centered on that element and
%   incorporates (window-1)/2 elements before and after the element.
%
%   Function is broken into two parts.  The 1d-2d solution, and the
%   n-dimensional solution.  The 1d-2d solution is the fastest that I have 
%   been able to come up with, whereas the n-dimensional solution trades
%   speed for versatility.
%
%   Function includes some code at the end so that the user can do their
%   own speed testing using the TIMEIT function.
%
%   Has been heavily tested in 1d-2d case, and lightly tested in 3d case.
%   Should work in n-dimensional space, but has not been tested in more 
%   than 3 dimensions other than to make sure it does not return an error.
%   Has not been tested for complex inputs.

%rounds even window sizes down to next lowest odd number
if mod(window,2)==0
    window=window-1;
end

%Calculates the number of elements in before and after the central element
%to incorporate in the moving mean.  Round command is just present to deal
%with the potential problem of division leaving a very small decimal, ie.
%2.000000000001.
halfspace=round((window-1)/2);

%calculates the size of the input data set
n=size(data);

%Computes the beginning and ending column for each moving
%average compuation.  Divide is the number of elements that are
%incorporated in each moving average.
start=[ones(1,halfspace+1) 2:(n(dim)-halfspace)];
stop=[(1+halfspace):n(dim) ones(1,halfspace)*n(dim)];
divide=stop-start+1;

%Calculates the moving average by calculating the sum of elements
%from the start row to the stop row for each central element,
%and then dividing by the number of elements used in that sum
%to get the average for that central element.
%Implemented by calculating the moving sum of the full data
%set.  Cumulative sum for each central element is calculated by
%subtracting the cumulative sum for the row before the start row
%from the cumulative sum for the stop row.  Row references are
%adusted to take into account the fact that you can now
%reference a row<1.  Divides the series of cumulative sums for
%by the number of elements in each sum to get the moving
%average.
CumulativeSum=cumsum(data);
temp_sum=CumulativeSum(stop,:)-CumulativeSum(max(start-1,1),:);
temp_sum((start==1),:)=bsxfun(@plus,temp_sum((start==1),:),data(1,:));
result=bsxfun(@rdivide,temp_sum,divide');

end