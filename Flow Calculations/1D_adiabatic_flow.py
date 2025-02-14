import math

def adiabatic_relations(R=None, gamma=1.4, q=None, M_upstream=None, M_downstream=None, T_upstream=None, T_downstream=None, \
        P_upstream=None, P_downstream=None, rho_upstream=None, rho_downstream=None, T0_upstream=None, T0_downstream=None, \
        P0_upstream=None, P0_downstream=None, rho0_upstream=None, rho0_downstream=None, T_star_upstream=None, T_star_downstream=None, \
        P_star_upstream=None, P_star_downstream=None, rho_star_upstream=None, rho_star_downstream=None, a_upstream=None, \
        a_downstream=None, a0_upstream=None, a0_downstream=None, a_star_upstream=None, a_star_downstream=None):

    """
    Computes adiabatic flow properties given various inputs.
    Computes flow properties for upstream and downstream.
    """

    Cp = (gamma * R)/(gamma - 1)
    adiabatic_inputs = [R, gamma, q, M_upstream, M_downstream, T_upstream, T_downstream, \
        P_upstream, P_downstream, rho_upstream, rho_downstream, T0_upstream, T0_downstream, \
        P0_upstream, P0_downstream, rho0_upstream, rho0_downstream, T_star_upstream, T_star_downstream, \
        P_star_upstream, P_star_downstream, rho_star_upstream, rho_star_downstream, a_upstream, \
        a_downstream, a0_upstream, a0_downstream, a_star_upstream, a_star_downstream]
    
    while None in adiabatic_inputs:
        
        # Compute Upstream Mach
        
        # Compute Downstream Mach
        
        # Compute missing upstream values
        if a_upstream is None and T_upstream is not None:
            a_upstream = math.sqrt(gamma * R * T_upstream)
        if a0_upstream is not None and T_star_upstream is not None and T0_upstream is not None:
            a_star_upstream = (a0_upstream * math.sqrt(T_star_upstream))/(math.sqrt(T0_upstream))
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
        
        # Compute missing downstream values
        if a_downstream is None and T_downstream is not None:
            a_downstream = math.sqrt(gamma * R * T_downstream)
        if a0_downstream is not None and T_star_downstream is not None and T0_downstream is not None:
            a_star_downstream = (a0_downstream * math.sqrt(T_star_downstream))/(math.sqrt(T0_downstream))
        if T0_downstream is not None:
            a0_downstream = math.sqrt(gamma * R * T0_downstream)
            
        if rho_downstream is None and P_upstream is not None and P_downstream is not None \
            and T_upstream is not None and T_downstream is not None and rho_upstream is not None:
            rho_downstream = 
        if P_downstream is None and rho_upstream is not None and rho_downstream is not None \
            and P_upstream is not None and T_upstream is not None and T_downstream is not None:
            P_downstream = 
        if T_downstream is None and rho_upstream is not None and rho_downstream is not None \
            and P_upstream is not None and P_downstream is not None and T_upstream is not None:
            T_downstream = 
        
        if P_upstream is None and P_downstream is not None and M_upstream is not None and M_downstream is not None:
            P_upstream = P_downstream / ((1 + gamma * (M_upstream ** 2)) / (1 + gamma * (M_downstream ** 2)))
        if T_upstream is None and T_downstream is not None and M_upstream is not None and M_downstream is not None:
            T_upstream = T_downstream / (((1 + gamma * (M_upstream ** 2)) / (1 + gamma * (M_downstream ** 2))) * \
                ((M_downstream / M_upstream) ** 2))
        if rho_upstream is None and rho_downstream is not None and M_upstream is not None and M_downstream is not None:
            rho_upstream = rho_downstream / (((1 + gamma * (M_downstream ** 2)) / (1 + gamma * (M_upstream ** 2))) * \
                ((M_upstream / M_downstream) ** 2))
        
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
            
            
            
        if M is not None:
            T_ratio = 1 + (gamma - 1) / 2 * M ** 2
            P_ratio = T_ratio ** (gamma / (gamma - 1))
            rho_ratio = T_ratio ** (1 / (gamma - 1))
            M_star = M / math.sqrt(T_ratio)
            T_star_ratio = 1 / T_ratio ** 2
            P_star_ratio = T_star_ratio ** (gamma / (gamma - 1))
            rho_star_ratio = T_star_ratio ** (1 / (gamma - 1))
        if T is None and T0 is not None and M is not None:
            T = T0 / T_ratio
        if P is None and P0 is not None and M is not None:
            P = P0 / P_ratio
        if rho is None and rho0 is not None and M is not None:
            rho = rho0 / rho_ratio
        if T0 is None and T is not None and M is not None:
            T0 = T * T_ratio
        if P0 is None and P is not None and M is not None:
            P0 = P * P_ratio
        if rho0 is None and rho is not None and M is not None:
            rho0 = rho * rho_ratio
        if M is not None and T0 is not None:
            T_star = T0 * T_star_ratio
        if M is not None and P0 is not None:
            P_star = P0 * P_star_ratio
        if M is not None and rho0 is not None:
            rho_star = rho0 * rho_star_ratio
        if P is None and rho is not None and P0 is not None:
            P = P0 * (rho / rho0) ** gamma
        if P is None and T is not None and P0 is not None:
            P = P0 * T * ((T / T0) ** (1 / (gamma - 1))) / T0
        if rho is None and T is not None and rho0 is not None:
            rho = rho0 * ((T * (T / T0) ** (1 / (gamma - 1))) / \
                T0) ** (1 / gamma)
        if rho is None and P is not None and rho0 is not None:
            rho = rho0 * (P / P0) ** (1 / gamma)
        if T is None and P is not None and T0 is not None:
            T = T0 / ((P0 / P) ** ((gamma - 1)/ gamma))
        if T is None and rho is not None and T0 is not None:
            T = T0 / (((rho0 / rho) ** gamma) ** ((gamma - 1) / gamma))
        
        adiabatic_inputs = [R, gamma, q, M_upstream, M_downstream, T_upstream, T_downstream, \
        P_upstream, P_downstream, rho_upstream, rho_downstream, T0_upstream, T0_downstream, \
        P0_upstream, P0_downstream, rho0_upstream, rho0_downstream, T_star_upstream, T_star_downstream, \
        P_star_upstream, P_star_downstream, rho_star_upstream, rho_star_downstream, a_upstream, \
        a_downstream, a0_upstream, a0_downstream, a_star_upstream, a_star_downstream]

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
        "Upstream Star Temperature Ratio (T*/T0)": T_star_ratio_upstream,
        "Downstream Star Temperature Ratio (T*/T0)": T_star_ratio_downstream,
        "Upstream Star Pressure Ratio (P*/P0)": P_star_ratio_upstream,
        "Downstream Star Pressure Ratio (P*/P0)": P_star_ratio_downstream,
        "Upstream Star Density Ratio (rho*/rho0)": rho_star_ratio_upstream,
        "Downstream Star Density Ratio (rho*/rho0)": rho_star_ratio_downstream,
        "Upstream Star Mach (M*)": M_star_upstream,
        "Downstream Star Mach (M*)": M_star_downstream,
        "Upstream Star Temperature (T*)": T_star_upstream,
        "Downstream Star Temperature (T*)": T_star_downstream,
        "Upstream Star Pressure (P*)": P_star_upstream,
        "Downstream Star Pressure (P*)": P_star_downstream,
        "Upstream Star Density (rho*)": rho_star_upstream,
        "Downstream Star Density (rho*)": rho_star_downstream,
        "Upstream Star Speed of Sound (a*)": a_star_upstream,
        "Downstream Star Speed of Sound (a*)": a_star_downstream
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
    P_star_upstream = input("Enter upstream star pressure (P*) or press enter to skip: ")
    T_star_upstream = input("Enter upstream star temperature (T*) or press enter to skip: ")
    rho_star_upstream = input("Enter upstream star density (rho*) or press enter: ")
    
    P_star_downstream = input("Enter downstream star pressure (P*) or press enter to skip: ")
    T_star_downstream = input("Enter downstream star temperature (T*) or press enter to skip: ")
    rho_star_downstream = input("Enter downstream star density (rho*) or press enter: ")
    
    a_star_upstream = input("Enter upstream star speed of sound (a*) or press enter: ")
    a_star_downstream = input("Enter downstream star speed of sound (a*) or press enter: ")
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

    P_star_upstream = float(P_star_upstream) if P_star_upstream else None
    T_star_upstream = float(T_star_upstream) if T_star_upstream else None
    rho_star_upstream = float(rho_star_upstream) if rho_star_upstream else None
    
    P_star_downstream = float(P_star_downstream) if P_star_downstream else None
    T_star_downstream = float(T_star_downstream) if T_star_downstream else None
    rho_star_downstream = float(rho_star_downstream) if rho_star_downstream else None
    
    a_upstream = float(a_upstream) if a_upstream else None
    a_downstream = float(a_downstream) if a_downstream else None
    a0_upstream = float(a0_upstream) if a0_upstream else None
    a0_downstream = float(a0_downstream) if a0_downstream else None
    a_star_upstream = float(a_star_upstream) if a_star_upstream else None
    a_star_downstream = float(a_star_downstream) if a_star_downstream else None

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
        P0_upstream, P0_downstream, rho0_upstream, rho0_downstream, T_star_upstream, T_star_downstream, \
        P_star_upstream, P_star_downstream, rho_star_upstream, rho_star_downstream, a_upstream, \
        a_downstream, a0_upstream, a0_downstream, a_star_upstream, a_star_downstream)

    for key, value in results.items():
        if value is not None:
            print(f"{key}: {value:.4f}")

if __name__ == "__main__":
    main()
