# MAC0210 - EP1
# Alunos:   Arthur Coser Marinho                NUSP: 
#           Ludmila Ferreira Vicente e Silva    NUSP: 7557136

function root = newton (f, df, x0)
    
    i = 0;
    err = 0;
    lim_err = 10^-3

    do
        prev = x0;
    
        x0 = x0 - (f (x0) ./ df (x0));
        

        err = abs (abs (prev - abs (x0)));

        i++;
    until (err < lim_err)
    
    root = x0;

endfunction
