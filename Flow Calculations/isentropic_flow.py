import math

def isentropic_relations(R=None, M=None, T=None, P=None, rho=None, a=None, gamma=1.4, \
    P0=None, T0=None, rho0=None, a0=None, P_star=None, T_star=None, rho_star=None, a_star=None):

    """
    Computes isentropic flow properties given various inputs.
    Optionally calculates absolute values if stagnation conditions are provided.
    """

    isentropic_inputs = [R, M, T, P, rho, a, gamma, P0, T0, rho0, \
        a0, P_star, T_star, rho_star, a_star]
    
    while None in isentropic_inputs:
        # Compute missing values if possible
        if T is None and P is not None and rho is not None:
            T = P / (rho * R)
        if rho is None and P is not None and T is not None:
            rho = P / (R * T)
        if P is None and rho is not None and T is not None:
            P = rho * R * T
        if a is None and T is not None:
            a = math.sqrt(gamma * R * T)
        if M is None and T is not None and T0 is not None:
            M = (math.sqrt(2 * (T0 - T))/(math.sqrt(T) * math.sqrt(gamma - 1)))
        if M is None and P is not None and P0 is not None:
            M = (math.sqrt(2) * math.sqrt(P0 * (P0 / P) ** (-1/gamma) - \
                P))/math.sqrt(gamma * P - P)
        if M is None and rho is not None and rho0 is not None:
            M = (math.sqrt(2) * math.sqrt(rho * (rho0 / rho) ** (gamma) - \
                rho0))/math.sqrt(gamma * rho0 - rho0)

        # Compute ratios (stagnation to actual values) from M
        if M is not None:
            T_ratio = 1 + (gamma - 1) / 2 * M ** 2
            P_ratio = T_ratio ** (gamma / (gamma - 1))
            rho_ratio = T_ratio ** (1 / (gamma - 1))

            # Compute star conditions from M
            M_star = M / math.sqrt(T_ratio)
            T_star_ratio = 1 / T_ratio ** 2
            P_star_ratio = T_star_ratio ** (gamma / (gamma - 1))
            rho_star_ratio = T_star_ratio ** (1 / (gamma - 1))

        # Compute absolute values if stagnation conditions and M are given
        if T is None and T0 is not None and M is not None:
            T = T0 / T_ratio
        if P is None and P0 is not None and M is not None:
            P = P0 / P_ratio
        if rho is None and rho0 is not None and M is not None:
            rho = rho0 / rho_ratio

        # Compute stagnation values if actual values and M are given
        if T0 is None and T is not None and M is not None:
            T0 = T * T_ratio
        if P0 is None and P is not None and M is not None:
            P0 = P * P_ratio
        if rho0 is None and rho is not None and M is not None:
            rho0 = rho * rho_ratio
        if T0 is not None:
            a0 = math.sqrt(gamma * R * T0)

        # Compute star values if stagnation values and M are given
        if M is not None and T0 is not None:
            T_star = T0 * T_star_ratio
        if M is not None and P0 is not None:
            P_star = P0 * P_star_ratio
        if M is not None and rho0 is not None:
            rho_star = rho0 * rho_star_ratio
        if a0 is not None and T_star is not None and T0 is not None:
            a_star = (a0 * math.sqrt(T_star))/(math.sqrt(T0))

        # Compute actual values from each other and M
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
        
        isentropic_inputs = [R, M, T, P, rho, a, gamma, P0, T0, rho0, \
        a0, P_star, T_star, rho_star, a_star]
        print(isentropic_inputs)

    return {
        "R": R,
        "Mach": M,
        "Temperature (T)": T,
        "Pressure (P)": P,
        "Density (rho)": rho,
        "Speed of Sound (a)": a,
        "Stagnation Temperature Ratio (T0/T)": T_ratio,
        "Stagnation Pressure Ratio (P0/P)": P_ratio,
        "Stagnation Density Ratio (rho0/rho)": rho_ratio,
        "Stagnation Temperature (T0)": T0,
        "Stagnation Pressure (P0)": P0,
        "Stagnation Density (rho0)": rho0,
        "Stagnation Speed of Sound (a0)": a0,
        "Star Temperature Ratio (T*/T0)": T_star_ratio,
        "Star Pressure Ratio (P*/P0)": P_star_ratio,
        "Star Density Ratio (rho*/rho0)": rho_star_ratio,
        "Star Mach (M*)": M_star,
        "Star Temperature (T*)": T_star,
        "Star Pressure (P*)": P_star,
        "Star Density (rho*)": rho_star,
        "Star Speed of Sound (a*)": a_star
    }

def main():
    R_unit = input("Enter units (SI or US): ")

    gamma = float(input("Enter specific heat ratio (gamma, default 1.4):") or 1.4)

    M = input("Enter Mach number (or press enter to skip): ")
    T = input("Enter temperature (or press enter to skip): ")
    P = input("Enter pressure (or press enter to skip): ")
    rho = input("Enter density (or press enter to skip): ")
    a = input("Enter speed of sound (or press enter to skip): ")

    # Optional stagnation conditions
    P0 = input("Enter stagnation pressure (P0) or press enter to skip: ")
    T0 = input("Enter stagnation temperature (T0) or press enter to skip: ")
    rho0 = input("Enter stagnation density (rho0) or press enter to skip: ")
    a0 = input("Enter stagnation speed of sound (a0) or press enter to skip: ")

    # Optional star conditions
    P_star = input("Enter star pressure (P*) or press enter to skip: ")
    T_star = input("Enter star temperature (T*) or press enter to skip: ")
    rho_star = input("Enter star density (rho*) or press enter: ")
    a_star = input("Enter star speed of sound (a*) or press enter: ")

    M = float(M) if M else None
    T = float(T) if T else None
    P = float(P) if P else None
    rho = float(rho) if rho else None
    a = float(a) if a else None

    P0 = float(P0) if P0 else None
    T0 = float(T0) if T0 else None
    rho0 = float(rho0) if rho0 else None
    a0 = float(a0) if a0 else None

    P_star = float(P_star) if P_star else None
    T_star = float(T_star) if T_star else None
    rho_star = float(rho_star) if rho_star else None
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

    results = isentropic_relations(R, M, T, P, rho, a, gamma, P0, T0, rho0, \
        a0, P_star, T_star, rho_star, a_star)

    for key, value in results.items():
        if value is not None:
            print(f"{key}: {value:.4f}")

if __name__ == "__main__":
    main()
