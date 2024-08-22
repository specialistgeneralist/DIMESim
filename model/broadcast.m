function output = broadcast(p)

Z = rand();

if Z <= p
    result = 1;
else
    result = -1;
end

output = result;
end