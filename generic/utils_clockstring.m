% convenience function for reading the clock
function s = utils_clockstring()
t = clock;
s = sprintf('%04d_%02d_%02d_%02d_%02d', t(1), t(2),t(3), t(4), t(5));
end
