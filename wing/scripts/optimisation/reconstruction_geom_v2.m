function x = reconstruction_geom_v2(x0)

    x.BoxGeo            = [0  x0(1)  x0(2);
                           1  x0(3)  x0(4)]; 
    x.tSkin             = [0   x0(5);
                           1   x0(6)];
    x.tWeb              = [0   x0(7);
                           1   x0(8)];
    x.StringerHeight    = [0   x0(10);
                           1   x0(11)];
    x.StringerThickness = [0   x0(12);
                           1   x0(13)];
    
    x.Stringer = x0(9);

end
