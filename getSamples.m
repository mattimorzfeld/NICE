function out = getSamples(nos,u,l)
out = u*(sqrt(l).*randn(length(l),nos));