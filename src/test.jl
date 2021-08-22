module A
    module B
        fact(n) = n == 0 ? 1 : fact(n-1) * n
        export fact
    end
    #using .B
end

        