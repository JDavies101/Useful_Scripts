import math

def adiabatic_relations(R=None, gamma=1.4, q=None, M_upstream=None, M_downstream=None, T_upstream=None, T_downstream=None, \
        P_upstream=None, P_downstream=None, rho_upstream=None, rho_downstream=None, T0_upstream=None, T0_downstream=None, \
        P0_upstream=None, P0_downstream=None, rho0_upstream=None, rho0_downstream=None, T_star=None, P_star=None, \
        rho_star=None, a_upstream=None, a_downstream=None, a0_upstream=None, a0_downstream=None, a_star=None):

    """
    Computes adiabatic flow properties given various inputs.
    Computes flow properties for upstream and downstream.
    """

    Cp = (gamma * R)/(gamma - 1)

    M_star = None
    T0_star = None
    P0_star = None
    P_ratio = None
    T_ratio = None  

    adiabatic_inputs = [R, gamma, q, M_upstream, M_downstream, T_upstream, T_downstream, \
        P_upstream, P_downstream, rho_upstream, rho_downstream, T0_upstream, T0_downstream, \
        P0_upstream, P0_downstream, rho0_upstream, rho0_downstream, T_star, P_star, rho_star, \
        a_upstream, a_downstream, a0_upstream, a0_downstream, a_star, M_star, T0_star, P0_star,\
        P_ratio, T_ratio]
    
    while None in adiabatic_inputs:
        
        # Compute Upstream Mach
        if M_upstream is None and M_downstream is not None and P_upstream is not None and P_downstream is not None:
            M_upstream = math.sqrt(gamma * (M_downstream ** 2) * P_downstream + P_downstream - P_upstream) / \
                (math.sqrt(gamma) * math.sqrt(P_upstream))
        if M_upstream is None and M_downstream is not None and T_upstream is not None and T_downstream is not None:
            M_upstream = math.sqrt((math.sqrt(T_downstream * (gamma * M_downstream ** 2 + 1) ** 2 * (-4 * T_upstream * gamma * M_downstream ** 2\
                        + gamma ** 2 * M_downstream ** 4 * T_downstream + 2 * gamma * M_downstream ** 2 * T_downstream + T_downstream)\
                        - 2 * T_upstream * gamma * M_downstream ** 2 + gamma ** 2 * M_downstream ** 4 * T_downstream + 2 * gamma \
                        * M_downstream ** 2 * T_downstream + T_downstream) / (T_upstream * gamma ** 2 * M_downstream ** 2)) / math.sqrt(2))
        if M_upstream is None and rho_upstream is not None and rho_downstream is not None:
            M_upstream = (M_downstream * math.sqrt(rho_downstream)) / math.sqrt(rho_upstream * gamma * M_downstream ** 2 + rho_upstream \
                        - gamma * M_downstream ** 2 * rho_downstream)
        
        if M_upstream is None and P_upstream is not None and P_star is not None:
            M_upstream = (math.sqrt(P_upstream - P_star * (gamma + 1)) * math.sqrt(-1 / gamma)) / math.sqrt(P_upstream)
        if M_upstream is None and rho_upstream is not None and rho_star is not None:
            M_upstream = math.sqrt(rho_star) / math.sqrt(rho_upstream * (gamma + 1) - rho_upstream * gamma)
        if M_upstream is None and T_upstream is not None and T_star is not None:
            M_upstream = math.sqrt((gamma * math.sqrt(T_star) * math.sqrt(gamma ** 2 * T_star - 4 * gamma * T_upstream + 2 * gamma * T_star + T_star) \
                        - math.sqrt(T_star) * math.sqrt(gamma ** 2 * T_star - 4 * gamma * T_upstream + 2 * gamma * T_star + T_star) + gamma ** 2 \
                        * T_star - 2 * gamma * T_upstream + 2 * gamma * T_star + T_star) / (gamma ** 2 * T_upstream)) / math.sqrt(2)
        
        if M_upstream is None and T_upstream is not None and T0_upstream is not None:
            M_upstream = math.sqrt((T0_upstream / T_upstream - 1) * 2 / (gamma - 1))

        if M_upstream is None and M_downstream is not None and T0_downstream is not None and T0_upstream is not None:
            M_upstream = math.sqrt(math.sqrt(T0_downstream) * math.sqrt(T0_downstream * (M_downstream ** 2 * gamma + 1) ** 2 - M_downstream ** 2 * T0_upstream\
                        * (gamma + 1) * (M_downstream ** 2 * gamma - M_downstream ** 2 + 2)) * (M_downstream ** 2 * gamma + 1) + T0_downstream \
                        * (M_downstream ** 2 * gamma + 1) ** 2 - M_downstream ** 2 * T0_upstream * gamma * (M_downstream ** 2 * gamma - M_downstream ** 2\
                        + 2)) * math.sqrt(-1 / (T0_downstream * (gamma - 1) * (M_downstream ** 2 * gamma - 1) ** 2 - M_upstream ** 2 * T0_downstream \
                        * gamma ** 2 *(M_downstream ** 2 * gamma - M_downstream ** 2 + 2)))

        if M_star is None and M_upstream is not None and T0_upstream is not None and T_upstream is not None:
            M_star = M_upstream / math.sqrt(T0_upstream / T_upstream)
        if M_upstream is None and M_star is not None and T0_upstream is not None and T_upstream is not None:
            M_upstream = M_star * math.sqrt(T0_upstream / T_upstream)

        # Compute Downstream Mach
        if M_downstream is None and M_upstream is not None and P_upstream is not None and P_downstream is not None:
            M_downstream = (math.sqrt(-1 / gamma) * math.sqrt(P_downstream - P_upstream * (gamma * (M_upstream\
                 ** 2) + 1))) / math.sqrt(P_downstream)
        if M_downstream is None and M_upstream is not None and T_upstream is not None and T_downstream is not None:
            M_downstream = math.sqrt((T_upstream * gamma ** 2 * M_upstream ** 4 - math.sqrt(T_upstream) * gamma * M_upstream ** 2\
                * math.sqrt(T_upstream * gamma ** 2 * M_upstream ** 4 + 2 * T_upstream * gamma * M_upstream ** 2 + T_upstream \
                - 4 * gamma * M_upstream ** 2 * T_downstream) - math.sqrt(T_upstream) * math.sqrt(T_upstream * gamma ** 2 * M_upstream ** 4 \
                + 2 * T_upstream * gamma * M_upstream ** 2 + T_upstream - 4 * gamma * M_upstream ** 2 * T_downstream) + 2 * T_upstream \
                * gamma * M_upstream ** 2 + T_upstream - 2 * gamma * M_upstream ** 2 * T_downstream) / (gamma ** 2 * M_upstream ** 2 \
                * T_downstream)) / math.sqrt(2)
        if M_downstream is None and rho_upstream is not None and rho_downstream is not None:
            M_downstream = (math.sqrt(rho_upstream) * M_upstream) / math.sqrt(-rho_upstream * gamma * M_upstream ** 2 + gamma * M_upstream ** 2 \
                            * rho_downstream + rho_downstream)
        
        if M_downstream is None and P_downstream is not None and P_star is not None:
            M_downstream = (math.sqrt(P_downstream - P_star * (gamma + 1)) * math.sqrt(-1 / gamma)) / math.sqrt(P_downstream)
        if M_downstream is None and rho_downstream is not None and rho_star is not None:
            M_downstream = math.sqrt(rho_star) / math.sqrt(rho_downstream * (gamma + 1) - rho_downstream * gamma)
        if M_downstream is None and T_downstream is not None and T_star is not None:
            M_downstream = math.sqrt((gamma * math.sqrt(T_star) * math.sqrt(gamma ** 2 * T_star - 4 * gamma * T_downstream + 2 * gamma * T_star + T_star) \
                        - math.sqrt(T_star) * math.sqrt(gamma ** 2 * T_star - 4 * gamma * T_downstream + 2 * gamma * T_star + T_star) + gamma ** 2 \
                        * T_star - 2 * gamma * T_downstream + 2 * gamma * T_star + T_star) / (gamma ** 2 * T_downstream)) / math.sqrt(2)
        
        if M_downstream is None and T_downstream is not None and T0_downstream is not None:
            M_downstream = math.sqrt((T0_downstream / T_downstream - 1) * 2 / (gamma - 1))

        if M_downstream is None and M_upstream is not None and T0_downstream is not None and T0_upstream is not None:
            M_downstream = math.sqrt(-1 / math.sqrt(M_upstream ** 4 * (T0_downstream - T0_upstream) * gamma ** 2 * (gamma - 1) + 2 * M_upstream ** 2 \
                            * (T0_downstream * gamma - T0_upstream * (gamma - 1)) * gamma - T0_upstream * (gamma - 1))) * math.sqrt(math.sqrt(\
                            M_upstream ** 4 * (T0_downstream * (gamma - 1) * (gamma + 1) - T0_upstream * gamma ** 2) + 2 * M_upstream ** 2 * (T0_downstream\
                            * (gamma + 1) - T0_upstream * gamma) - T0_upstream) * math.abs(M_upstream ** 2 * gamma + 1) * math.sqrt(-T0_upstream) \
                            + M_upstream ** 4 * (T0_downstream * (gamma - 1) - T0_upstream * gamma) * gamma + 2 * M_upstream ** 2 * (T0_downstream\
                            - T0_upstream) * gamma - T0_upstream)

        if M_star is None and M_downstream is not None and T0_downstream is not None and T_downstream is not None:
            M_star = M_downstream / math.sqrt(T0_downstream / T_downstream)
        if M_downstream is None and M_star is not None and T0_downstream is not None and T_downstream is not None:
            M_downstream = M_star * math.sqrt(T0_downstream / T_downstream) 
        
        # Compute missing upstream values
        if a_upstream is None and T_upstream is not None:
            a_upstream = math.sqrt(gamma * R * T_upstream)
        if a0_upstream is not None and T_star is not None and T0_upstream is not None:
            a_star = (a0_upstream * math.sqrt(T_star))/(math.sqrt(T0_upstream))
        if T0_upstream is not None:
            a0_upstream = math.sqrt(gamma * R * T0_upstream)
            
        if rho_upstream is None and P_upstream is not None and P_downstream is not None \
            and T_upstream is not None and T_downstream is not None:
            rho_upstream = (P_downstream * T_upstream) / (rho_downstream * P_upstream * T_downstream)
        if P_upstream is None and rho_upstream is not None and rho_downstream is not None \
            and P_downstream is not None and T_upstream is not None and T_downstream is not None:
            P_upstream = (rho_upstream * P_downstream * T_upstream) / (rho_downstream * T_downstream)
        if T_upstream is None and rho_upstream is not None and rho_downstream is not None \
            and P_upstream is not None and P_downstream is not None and T_downstream is not None:
            T_upstream = (rho_downstream * P_upstream * T_downstream) / (rho_upstream * P_downstream)
        
        if P_upstream is None and P_downstream is not None and M_upstream is not None and M_downstream is not None:
            P_upstream = P_downstream / ((1 + gamma * (M_upstream ** 2)) / (1 + gamma * (M_downstream ** 2)))
        if T_upstream is None and T_downstream is not None and M_upstream is not None and M_downstream is not None:
            T_upstream = T_downstream / (((1 + gamma * (M_upstream ** 2)) / (1 + gamma * (M_downstream ** 2))) * \
                ((M_downstream / M_upstream) ** 2))
        if rho_upstream is None and rho_downstream is not None and M_upstream is not None and M_downstream is not None:
            rho_upstream = rho_downstream / (((1 + gamma * (M_downstream ** 2)) / (1 + gamma * (M_upstream ** 2))) * \
                ((M_upstream / M_downstream) ** 2))
        
        if T_upstream is None and T0_upstream is not None and M_upstream is not None:
            T_upstream = T0_upstream / (1 + (gamma - 1) / 2 * M_upstream ** 2)
            
        if T0_upstream is None and T0_downstream is not None and q is not None:
            T0_upstream = T0_downstream - q / Cp
        if T0_upstream is None and T_upstream is not None and M_upstream is not None:
            T0_upstream = T_upstream * (1 + (gamma - 1) / 2 * M_upstream ** 2)
        
        if P0_upstream is None and T_upstream is not None and M_upstream is not None:
            P0_upstream = (T_upstream * (1 + (gamma - 1) / 2 * M_upstream ** 2)) ** (gamma / (gamma - 1))
        
        if rho0_upstream is None and T_upstream is not None and M_upstream is not None:
            rho0_upstream = (T_upstream * (1 + (gamma - 1) / 2 * M_upstream ** 2)) ** (1 / (gamma - 1))
        
        if T_star is None and T_upstream is not None and M_upstream is not None:
            T_star = T_upstream / ((M_upstream ** 2) * (1 + gamma) / 1 + gamma * (M_upstream ** 2))
        if P_star is None and P_upstream is not None and M_upstream is not None:
            P_star = P_upstream / ((1 + gamma) / 1 + gamma * (M_upstream ** 2))
        if rho_star is None and rho_upstream is not None and M_upstream is not None:
            rho_star = rho_upstream / M_upstream * ((1 + gamma * (M_upstream ** 2)) / (1 + gamma))
        
        if T0_star is None and T0_upstream is not None and M_upstream is not None:
            T0_star = T0_upstream / (((gamma + 1) * M_upstream ** 2) / (1 + gamma * (M_upstream ** 2)) * (2 + (gamma - 1) * (M_upstream ** 2)))
        if P0_star is None and P0_upstream is not None and M_upstream is not None:
            P0_star = P0_upstream / (((1 + gamma) / (1 + gamma * (M_upstream ** 2))) * ((2 + (gamma - 1) * M_upstream ** 2) / (gamma + 1)) ** (gamma / (gamma - 1)))
        
        # Compute missing downstream values
        if a_downstream is None and T_downstream is not None:
            a_downstream = math.sqrt(gamma * R * T_downstream)
        if a0_downstream is not None and T_star is not None and T0_downstream is not None:
            a_star = (a0_downstream * math.sqrt(T_star))/(math.sqrt(T0_downstream))
        if T0_downstream is not None:
            a0_downstream = math.sqrt(gamma * R * T0_downstream)
            
        if rho_downstream is None and P_upstream is not None and P_downstream is not None \
            and T_upstream is not None and T_downstream is not None and rho_upstream is not None:
            rho_downstream = (rho_upstream * P_downstream * T_upstream) / (P_upstream * T_downstream)
        if P_downstream is None and rho_upstream is not None and rho_downstream is not None \
            and P_upstream is not None and T_upstream is not None and T_downstream is not None:
            P_downstream = (rho_downstream * P_upstream * T_downstream) / (rho_upstream * T_upstream)
        if T_downstream is None and rho_upstream is not None and rho_downstream is not None \
            and P_upstream is not None and P_downstream is not None and T_upstream is not None:
            T_downstream = (rho_upstream * P_downstream * T_upstream) / (rho_downstream * P_upstream)
        
        if P_downstream is None and P_upstream is not None and M_upstream is not None and M_downstream is not None:
            P_downstream = P_upstream * ((1 + gamma * (M_upstream ** 2)) / (1 + gamma * (M_downstream ** 2)))
        if T_downstream is None and T_upstream is not None and M_upstream is not None and M_downstream is not None:
            T_downstream = T_upstream * (((1 + gamma * (M_upstream ** 2)) / (1 + gamma * (M_downstream ** 2))) * \
                ((M_downstream / M_upstream) ** 2))
        if rho_downstream is None and rho_upstream is not None and M_upstream is not None and M_downstream is not None:
            rho_downstream = rho_upstream * (((1 + gamma * (M_downstream ** 2)) / (1 + gamma * (M_upstream ** 2))) * \
                ((M_upstream / M_downstream) ** 2))
        
        if T_downstream is None and T0_downstream is not None and M_downstream is not None:
            T_downstream = T0_downstream / (1 + (gamma - 1) / 2 * M_downstream ** 2)
        
        if T0_downstream is None and T0_upstream is not None and q is not None:
            T0_downstream = q / Cp - T0_upstream
        if T0_downstream is None and T_downstream is not None and M_downstream is not None:
            T0_downstream = T_downstream * (1 + (gamma - 1) / 2 * M_downstream ** 2)
        
        if P0_downstream is None and T_downstream is not None and M_downstream is not None:
            P0_downstream = (T_downstream * (1 + (gamma - 1) / 2 * M_downstream ** 2)) ** (gamma / (gamma - 1))
        
        if rho0_downstream is None and T_downstream is not None and M_downstream is not None:
            rho0_downstream = (T_downstream * (1 + (gamma - 1) / 2 * M_downstream ** 2)) ** (1 / (gamma - 1))
        
        if T_star is None and T_downstream is not None and M_downstream is not None:
            T_star = T_downstream / ((M_downstream ** 2) * (1 + gamma) / 1 + gamma * (M_downstream ** 2))
        if P_star is None and P_downstream is not None and M_downstream is not None:
            P_star = P_downstream / ((1 + gamma) / 1 + gamma * (M_downstream ** 2))
        if rho_star is None and rho_downstream is not None and M_downstream is not None:
            rho_star = rho_downstream / M_downstream * ((1 + gamma * (M_downstream ** 2)) / (1 + gamma))
        
        if T0_star is None and T0_downstream is not None and M_downstream is not None:
            T0_star = T0_downstream / (((gamma + 1) * M_downstream ** 2) / (1 + gamma * (M_downstream ** 2)) \
                * (2 + (gamma - 1) * (M_downstream ** 2)))
        if P0_star is None and P0_downstream is not None and M_downstream is not None:
            P0_star = P0_downstream / (((1 + gamma) / (1 + gamma * (M_downstream ** 2))) * ((2 + (gamma - 1) \
                * M_downstream ** 2) / (gamma + 1)) ** (gamma / (gamma - 1)))
        
        # Compute missing flow values
        if q is None and T0_upstream is not None and T0_downstream is not None:
            q = Cp * (T0_downstream - T0_upstream)
        
        if P_ratio is None and M_upstream is not None and M_downstream is not None:
            P_ratio = (1 + gamma * (M_upstream ** 2)) / (1 + gamma * (M_downstream ** 2)) * \
                ((1 + (gamma - 1) / 2 * (M_downstream ** 2)) / (1 + (gamma - 1) * (M_upstream ** 2))) \
                ** (gamma / (gamma - 1))
        if T_ratio is None and M_upstream is not None and M_downstream is not None:
            T_ratio = (((1 + gamma * (M_upstream ** 2)) / (1 + gamma * (M_downstream ** 2))) ** 2) * \
                ((M_downstream / M_upstream) ** 2) * ((1 + (gamma - 1) / 2 * (M_downstream ** 2)) / \
                (1 + (gamma - 1) * (M_upstream ** 2)))
        
        adiabatic_inputs = [R, gamma, q, M_upstream, M_downstream, T_upstream, T_downstream, \
        P_upstream, P_downstream, rho_upstream, rho_downstream, T0_upstream, T0_downstream, \
        P0_upstream, P0_downstream, rho0_upstream, rho0_downstream, T_star, P_star, rho_star, \
        a_upstream, a_downstream, a0_upstream, a0_downstream, a_star, M_star, T0_star, P0_star,\
        P_ratio, T_ratio]

        print(adiabatic_inputs)

    return {
        "R": R,
        "Cp": Cp,
        "Upstream Mach": M_upstream,
        "Downstream Mach": M_downstream,
        "Upstream Temperature (T)": T_upstream,
        "Downstream Temperature (T)": T_downstream,
        "Upstream Pressure (P)": P_upstream,
        "Downstream Pressure (P)": P_downstream,
        "Upstream Density (rho)": rho_upstream,
        "Downstream Density (rho)": rho_downstream,
        "Upstream Speed of Sound (a)": a_upstream,
        "Downstream Speed of Sound (a)": a_downstream,
        "Stagnation Temperature Ratio (T02/T01)": T_ratio,
        "Stagnation Pressure Ratio (P02/P01)": P_ratio,
        "Upstream Stagnation Temperature (T0)": T0_upstream,
        "Downstream Stagnation Temperature (T0)": T0_downstream,
        "Upstream Stagnation Pressure (P0)": P0_upstream,
        "Downstream Stagnation Pressure (P0)": P0_downstream,
        "Upstream Stagnation Density (rho0)": rho0_upstream,
        "Downstream Stagnation Density (rho0)": rho0_downstream,
        "Upstream Stagnation Speed of Sound (a0)": a0_upstream,
        "Downstream Stagnation Speed of Sound (a0)": a0_downstream,
        "Star Mach (M*)": M_star,
        "Star Temperature (T*)": T_star,
        "Star Pressure (P*)": P_star,
        "Star Density (rho*)": rho_star,
        "Star Speed of Sound (a*)": a_star,
        "Stagnation Star Temperature (T0*)": T0_star,
        "Stagnation Star Pressure (P0*)": P0_star
    }

def main():
    R_unit = input("Enter units (SI or US): ")

    gamma = float(input("Enter specific heat ratio (gamma, default 1.4):") or 1.4)
    
    q = input("Enter heat added per unit mass (or press enter to skip): ")

    M_upstream = input("Enter upstream Mach number (or press enter to skip): ")
    T_upstream = input("Enter upstream temperature (or press enter to skip): ")
    P_upstream = input("Enter upstream pressure (or press enter to skip): ")
    rho_upstream = input("Enter upstream density (or press enter to skip): ")
    
    M_downstream = input("Enter downstream Mach number (or press enter to skip): ")
    T_downstream = input("Enter downstream temperature (or press enter to skip): ")
    P_downstream = input("Enter downstream pressure (or press enter to skip): ")
    rho_downstream = input("Enter downstream density (or press enter to skip): ")

    # Optional stagnation conditions
    P0_upstream = input("Enter upstream stagnation pressure (P0) or press enter to skip: ")
    T0_upstream = input("Enter upstream stagnation temperature (T0) or press enter to skip: ")
    rho0_upstream = input("Enter upstream stagnation density (rho0) or press enter to skip: ")
    
    P0_downstream = input("Enter downstream stagnation pressure (P0) or press enter to skip: ")
    T0_downstream = input("Enter downstream stagnation temperature (T0) or press enter to skip: ")
    rho0_downstream = input("Enter downstream stagnation density (rho0) or press enter to skip: ")
    
    # Optional star conditions
    P_star = input("Enter star pressure (P*) or press enter to skip: ")
    T_star = input("Enter star temperature (T*) or press enter to skip: ")
    rho_star = input("Enter star density (rho*) or press enter: ")
    
    a_star= input("Enter star speed of sound (a*) or press enter: ")
    a_upstream = input("Enter upstream speed of sound (or press enter to skip): ")
    a_downstream = input("Enter downstream speed of sound (or press enter to skip): ")
    a0_upstream = input("Enter upstream stagnation speed of sound (a0) or press enter to skip: ")
    a0_downstream = input("Enter downstream stagnation speed of sound (a0) or press enter to skip: ")

    # None assignment if no input given
    q = float(q) if q else None
    
    M_upstream = float(M_upstream) if M_upstream else None
    T_upstream = float(T_upstream) if T_upstream else None
    P_upstream = float(P_upstream) if P_upstream else None
    rho_upstream = float(rho_upstream) if rho_upstream else None
    
    M_downstream = float(M_downstream) if M_downstream else None
    T_downstream = float(T_downstream) if T_downstream else None
    P_downstream = float(P_downstream) if P_downstream else None
    rho_downstream = float(rho_downstream) if rho_downstream else None

    P0_upstream = float(P0_upstream) if P0_upstream else None
    T0_upstream = float(T0_upstream) if T0_upstream else None
    rho0_upstream = float(rho0_upstream) if rho0_upstream else None
    
    P0_downstream = float(P0_downstream) if P0_downstream else None
    T0_downstream = float(T0_downstream) if T0_downstream else None
    rho0_downstream = float(rho0_downstream) if rho0_downstream else None

    P_star = float(P_star) if P_star else None
    T_star = float(T_star) if T_star else None
    rho_star = float(rho_star) if rho_star else None
    
    a_upstream = float(a_upstream) if a_upstream else None
    a_downstream = float(a_downstream) if a_downstream else None
    a0_upstream = float(a0_upstream) if a0_upstream else None
    a0_downstream = float(a0_downstream) if a0_downstream else None
    a_star = float(a_star) if a_star else None

    # Determine specific gas constant
    if R_unit == "US":
        R = float(1716) # Specific gas constant for air in ft*lbf/(slug*R)
    elif R_unit == "SI":
        R = float(287.05)  # Specific gas constant for air in J/(kg*K)
    elif isinstance(R_unit, float):
        R = R_unit
    else:
        raise ValueError("Incorrect R input")

    results = adiabatic_relations(R, gamma, q, M_upstream, M_downstream, T_upstream, T_downstream, \
        P_upstream, P_downstream, rho_upstream, rho_downstream, T0_upstream, T0_downstream, \
        P0_upstream, P0_downstream, rho0_upstream, rho0_downstream, T_star, P_star, rho_star, \
        a_upstream, a_downstream, a0_upstream, a0_downstream, a_star)

    for key, value in results.items():
        if value is not None:
            print(f"{key}: {value:.4f}")

if __name__ == "__main__":
    main()
