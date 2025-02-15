from sympy import symbols, Eq, solve

R, gamma, q, M_upstream, M_downstream, T_upstream, T_downstream, P_upstream, P_downstream, rho_upstream, rho_downstream,\
        P0_upstream, P0_downstream, rho0_upstream, rho0_downstream, T_star, P_star, rho_star, a_upstream, a_downstream, \
        a0_upstream, a0_downstream, a_star, Cp, P_ratio, T_ratio, rho_ratio, P0_ratio, T0_ratio, P_star_ratio_upstream, \
        P_star_ratio_downstream, T_star_ratio_upstream, T_star_ratio_downstream, rho_star_ratio_upstream, rho_star_ratio_downstream,\
        P0_star_ratio_upstream, P0_star_ratio_downstream, T0_star_ratio_upstream, T0_downstream, T0_upstream\
            = symbols("R gamma q M_upstream M_downstream T_upstream T_downstream P_upstream P_downstream rho_upstream rho_downstream\
                      P0_upstream P0_downstream rho0_upstream rho0_downstream T_star P_star rho_star a_upstream a_downstream a0_upstream\
                      a0_downstream a_star Cp P_ratio T_ratio rho_ratio P0_ratio T0_ratio P_star_ratio_upstream P_star_ratio_downstream\
                      T_star_ratio_upstream T_star_ratio_downstream rho_star_ratio_upstream rho_star_ratio_downstream P0_star_ratio_upstream\
                      P0_star_ratio_downstream T0_star_ratio_upstream T0_downstream T0_upstream")

# Equation definitions
equations = [
    #q_eqn 
    Eq(Cp * (T0_downstream - T0_upstream), q),
    #P_ratio_eqn
    Eq((1 + gamma * M_upstream ** 2) / (1 + gamma * M_downstream ** 2), P_ratio),
    #T_ratio_eqn
    Eq(((1 + gamma * M_upstream ** 2) / (1 + gamma * M_downstream ** 2)) ** 2 * (M_downstream / M_upstream) ** 2, T_ratio),
    #rho_ratio_eqn
    Eq(((1 + gamma * M_downstream ** 2) / (1 + gamma * M_upstream ** 2)) ** 2 * (M_upstream / M_downstream) ** 2, rho_ratio),
    #P0_ratio_eqn
    Eq((1 + gamma * M_upstream ** 2) / (1 + gamma * M_downstream ** 2) * ((1 + (gamma - 1) / 2 * M_downstream ** 2) \
                        / (1 + (gamma - 1) / 2 * M_upstream ** 2)) ** (gamma / (gamma - 1)), P0_ratio),
    #T0_ratio_eqn
    Eq(((1 + gamma * M_upstream ** 2) / (1 + gamma * M_downstream ** 2)) ** 2 * (M_downstream / M_upstream) ** 2 \
                        * (1 + (gamma - 1) / 2 * M_downstream ** 2) / (1 + (gamma - 1) / 2 * M_upstream ** 2), T0_ratio),
    #P_star_ratio_upstream_eqn
    Eq((1 + gamma) / (1 + gamma * M_upstream ** 2), P_star_ratio_upstream),
    #P_star_ratio_downstream_eqn
    Eq((1 + gamma) / (1 + gamma * M_downstream ** 2), P_star_ratio_downstream),
    #T_star_ratio_upstream_eqn
    Eq(M_upstream ** 2 * ((1 + gamma) / (1 + gamma * M_upstream ** 2)) ** 2, T_star_ratio_upstream),
    #T_star_ratio_downstream_eqn
    Eq(M_downstream ** 2 * ((1 + gamma) / (1 + gamma * M_downstream ** 2)) ** 2, T_star_ratio_downstream),
    #rho_star_ratio_upstream_eqn
    Eq(1 / M_upstream ** 2 * (1 + gamma * M_upstream ** 2) / (1 + gamma), rho_star_ratio_upstream),
    #rho_star_ratio_downstream_eqn
    Eq(1 / M_downstream ** 2 * (1 + gamma * M_downstream ** 2) / (1 + gamma), rho_star_ratio_downstream),
    #P0_star_ratio_upstream_eqn
    Eq((1 + gamma) / (1+ gamma * M_upstream ** 2) * ((2 + (gamma - 1) * M_upstream ** 2) / (gamma + 1)) ** \
                                    (gamma / (gamma - 1)), P0_star_ratio_upstream),
    #P0_star_ratio_downstream_eqn
    Eq((1 + gamma) / (1+ gamma * M_downstream ** 2) * ((2 + (gamma - 1) * M_downstream ** 2) / (gamma + 1)) ** \
                                    (gamma / (gamma - 1)), P0_star_ratio_downstream),
    #T0_star_ratio_upstream_eqn
    Eq(((gamma + 1) * M_upstream ** 2) / (1 + gamma * M_upstream ** 2) ** 2 * (2 +(gamma - 1) * M_upstream ** 2), \
                                    T0_star_ratio_upstream),
    #T0_star_ratio_downstream_eqn
    Eq(((gamma + 1) * M_downstream ** 2) / (1 + gamma * M_downstream ** 2) ** 2 * (2 +(gamma - 1) * M_downstream \
                                    ** 2), T0_star_ratio_upstream),
    #ideal_upstream_eqn
    Eq(rho_upstream * R * T_upstream, P_upstream),
    #ideal_downstream_eqn
    Eq(rho_downstream * R * T_downstream, P_downstream),
    #Cp_eqn
    Eq(Cp, (gamma * R) / (gamma - 1)),
]   

def adiabatic_relations(eqs, known_values):

    """
    Computes adiabatic flow properties given various inputs.
    Computes flow properties for upstream and downstream.
    """
    all_variables = set().union(*[eq.free_symbols for eq in eqs])
    unknowns = list(all_variables - set(known_values.keys()))

    while unknowns:
        solved_count = 0 # Track number of unknowns solved in each iteration
        for eq in eqs:
            subbed_eqs = eq.subs(known_values)

            free_vars = list(subbed_eqs.free_symbols)
            # Try solving for unknowns
            if len(free_vars) == 1:
                solutions = solve(subbed_eqs, free_vars[0])

                # Filter out complex solutions
                real_solutions = [s.evalf() for s in solutions if s.is_real]
                positive_solutions = [s for s in real_solutions if s >= 0]

                # Filter out real and positive solutions
                if len(positive_solutions) == 1:
                    chosen_solution = positive_solutions[0]
                elif len(positive_solutions) > 1:
                    print(f"\nMultiple solutions found for {free_vars[0]}: {positive_solutions}")
                    chosen_solution = float(input(f"Enter preferred solution for {free_vars[0]}: "))
                elif real_solutions:  # If no positive, take any real solution
                    chosen_solution = real_solutions[0]
                else:
                    continue  # No valid solutions

                known_values[free_vars[0]] = chosen_solution
                unknowns.remove(free_vars[0])
                solved_count += 1
        
        if solved_count == 0:
            break  # Stop if no new unknowns were solved

    return known_values

def main():

    knowns = {}

    for var in [R, gamma, q, M_upstream, M_downstream, T_upstream, T_downstream, P_upstream, P_downstream, rho_upstream, rho_downstream,\
                P0_upstream, P0_downstream, rho0_upstream, rho0_downstream, T_star, P_star, rho_star, a_upstream, a_downstream, a0_upstream,\
                a0_downstream, a_star]:
        inputs = input(f"Enter {var} or press enter to skip: ")
        if inputs.strip():
            knowns[var] = float(inputs)
    
    print(knowns)
    
    results = adiabatic_relations(equations, knowns)

    print(results)

    if results:
        print("\nSolved Values:")
        for var, value in results.items():
            print(f"{var} = {value:.4f}")

        
if __name__ == "__main__":
    main()
