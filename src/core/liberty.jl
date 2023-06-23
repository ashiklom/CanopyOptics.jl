"""
liberty(D, xu, thick, baseline, element, c_factor, w_factor, l_factor, p_factor, optis)

Leaf optical properties (reflectance and transmittance) based on LIBERTY model.

Arguments are as follows:

----------------------------------------------------------------------
  Variable        Description                Typical values(and range) 
----------------------------------------------------------------------
  Cell Diameter(D)   Average leaf cell diameter (m-6)      40 (20 - 100)
----------------------------------------------------------------------
  Intercellular       Determinant for the amount of                        
  air space (xu)   radiative flux passing between cells 0.045 (0.01 - 0.1)
----------------------------------------------------------------------
  Leaf thickness(thick)  Arbitrary value to determine single                 
                          leaf reflectance and transmittance                              
                         from infinite reflectance criteria 1.6 (1 - 10)  
----------------------------------------------------------------------
  Baseline(baseline) Wavelength independent absorption to                
   absorption        compensate for changes in absolute   Fresh: 0.0006  
  	                     reflectance                    Dry: 0.0004  
----------------------------------------------------------------------
  Albino          Absorption in the visible region due                  
  absorption(element)      to lignin                      2 (0 - 4)    
----------------------------------------------------------------------
  Chlorophyll              Chlorophyll (pigment)                                          
  content (c_factor)         content (mg.m-2)             200 (0 - 600) 
----------------------------------------------------------------------
   Water                                                                
  content(w_factor)         Water content (g.m-2)         100 (0 - 500)
----------------------------------------------------------------------
  Lignin and                                                            
  Cellulose            Combined lignin and cellulose                         
  content (l_factor)        content (g.m-2)               40 (10 - 80)  
----------------------------------------------------------------------
  Nitrogen              Nitrogen content                                      
  content (p_factor)        (g.m-2)                       1 (0.3 - 2.0) 
----------------------------------------------------------------------
"""
function liberty(
        D::FT, xu::FT, thick::FT, baseline::FT,
        element::FT, c_factor::FT, w_factor::FT, l_factor::FT, p_factor::FT,
        optis
) where {FT <: Real}
    (k_chloro, k_water, ke, k_ligcell, k_prot) = optis
    coeff = @. D * (baseline + 
        k_chloro*c_factor +
        k_water*w_factor +
        ke*element +
        k_ligcell*l_factor +
        k_prot*p_factor
    )

    N = length(coeff)
    wavelength = zeros(N)
    trans = zeros(N)
    refl = zeros(N)

    for i in 1:N
        wavelength[i], trans[i], refl[i] = liberty_body(coeff[i], i, thick, xu)
    end

    return wavelength, trans, refl
end

function liberty_body(coeff, i, thick, xu)
    # Refractive index calculation
    N1 = 1.4891-(0.0005*(i-1))
    N0 = 1.0
    in_angle = 59
    # % /* Index of refraction */
    alpha = in_angle * π/180;
    beta = asin((N0/N1)*sin(alpha));
    # Values are computed but not used
    # para_r = (tan(alpha-beta))/(tan(alpha+beta));
    # vert_r = -(sin(alpha-beta))/(sin(alpha+beta));

    me = 0
    width = π/180
    for p in range(1, 90)
        alpha = p*π/180
        beta = asin(N0/N1*sin(alpha))
        plus = alpha+beta
        dif = alpha-beta
        refl = 0.5*((sin(dif))^2/(sin(plus))^2+(tan(dif))^2/(tan(plus))^2)
        me = me+(refl*sin(alpha)*cos(alpha)*width)
    end

    mi=0
    mint=0
    width=pi/180
    critical=asin(N0/N1)*180/pi
    for j=1:critical
        alpha=j*pi/180
        beta=asin(N0/N1*sin(alpha))
        plus=alpha+beta
        dif=alpha-beta
        refl=0.5*((sin(dif))^2/(sin(plus))^2+(tan(dif))^2/(tan(plus))^2)
        mint=mint+refl*sin(alpha)*cos(alpha)*width
    end
    mi=(1-sin(critical*pi/180)^2)+2*mint

    M=2*(1-(coeff+1)*exp(-coeff))/coeff^2
    T=((1-mi)*M)/(1-(mi*M))
    x=xu/(1-(1-2*xu)*T)
    a=me*T+x*T-me-T-x*me*T
    b=1+x*me*T-2*x^2*me^2*T
    c=2*me*x^2*T-x*T-2*x*me
    # Initial guess
    R=0.5
    for _ in 1:50
        next_R=-(a*R^2+c)/b
        R=next_R
    end

    # 	  /* The next bit works out transmittance based upon Benford...*/
    # 	  /* setting up unchanging parameters... */
    rb=(2*x*me)+(x*T)-(x*T*2*x*me)
    tb=sqrt(((R-rb)*(1-(R*rb)))/R)
    whole=trunc(thick)
    fraction=thick-whole
    # 	/* The next bit works out the fractional value */
    # 	/* for the interval between 1 and 2... */
    top=tb^(1+fraction)*((((1+tb)^2)-rb^2)^(1-fraction))
    bot1=(1+tb)^(2*(1-fraction))-rb^2
    bot2=1+(64/3)*fraction*(fraction-0.5)*(fraction-1)*0.001
    tif=top/(bot1*bot2)
    rif=(1+rb^2-tb^2-sqrt((1+rb^2-tb^2)^2-4*rb^2*(1-tif^2)))/(2*rb);    
    #   /* Now to work out for integral thickness greater than 2 ... */
    if (whole>=2) 
        prev_t=1
        prev_r=0
        for _ in 1:(whole-1)
            cur_t=(prev_t*tb)/(1-(prev_r*rb))
            cur_r=prev_r+(((prev_t*prev_t)*rb)/(1-(prev_r*rb)))
            prev_t=cur_t
            prev_r=cur_r
        end
    else
        cur_t=1
        cur_r=0
    end
    trans=(cur_t*tif)/(1-(rif*cur_r))
    refl=cur_r+((cur_t^2*rif)/(1-(rif*cur_r)))
    wavelength = 400 + ((i-1)*5)

    return wavelength, trans, refl
end
