# MAC0210 - EP1
# Alunos:   Arthur Coser Marinho                NUSP: 7210629 
#           Ludmila Ferreira Vicente e Silva    NUSP: 7557136

function b = newton_basins (f, df, l, u, p)
# mudei o cabecalho da funcao pq tava liberado no paca
    
    range_x = linspace (-l, l, p);
    range_y = linspace (-u, u, p);  
    [X, Y] = meshgrid (range_x, range_y);
    Z = X + Y * i;
    Z = newton (f, df, Z)
    
endfunction
