# MAC0210 - EP1
# Alunos:   Arthur Coser Marinho                NUSP: 
#           Ludmila Ferreira Vicente e Silva    NUSP: 7557136

function x0 = newton (f, df, x0)
    
    i = 0;
    err = 0;
    lim_err = 10^-3

    do
        prev = x0;
    
        x0 = x0 - (f (x0) ./ df (x0));
        
        i++;

        err = abs (abs (prev - abs (x0)));

    until (err < lim_err)
    
    disp(x0);

endfunction
