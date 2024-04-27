function X = fftfiltfilt(B,X)

w = length(B);
t = size(X,1);
X = [repmat(2*X(1,:),w,1)-X((w+1):-1:2,:); X; repmat(2*X(t,:),w,1)-X((t-1):-1:t-w,:)];

X = fftfilt(B,X);
X = X(end:-1:1,:);

X = fftfilt(B,X);
X = X(end:-1:1,:);

X([1:w t+w+(1:w)],:) = [];

end