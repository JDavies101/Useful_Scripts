from sympy import symbols, Eq, solve, nsolve, sqrt

R, gamma, q, M_upstream, M_downstream, T_upstream, T_downstream, P_upstream, P_downstream, rho_upstream, rho_downstream,\
    P0_upstream, P0_downstream, rho0_upstream, rho0_downstream, T_star, P_star, rho_star, a_upstream, a_downstream, \
    a0_upstream, a0_downstream, a_star, Cp, P_ratio, T_ratio, rho_ratio, P0_ratio, T0_ratio, P_star_ratio_upstream, \
    P_star_ratio_downstream, T_star_ratio_upstream, T_star_ratio_downstream, rho_star_ratio_upstream, rho_star_ratio_downstream,\
    P0_star_ratio_upstream, P0_star_ratio_downstream, T0_star_ratio_upstream, T0_star_ratio_downstream, T0_downstream, T0_upstream, P0_star, T0_star\
        = symbols("R, gamma, q, M_upstream, M_downstream, T_upstream, T_downstream, P_upstream, P_downstream, rho_upstream, rho_downstream,\
    P0_upstream, P0_downstream, rho0_upstream, rho0_downstream, T_star, P_star, rho_star, a_upstream, a_downstream, \
    a0_upstream, a0_downstream, a_star, Cp, P_ratio, T_ratio, rho_ratio, P0_ratio, T0_ratio, P_star_ratio_upstream, \
    P_star_ratio_downstream, T_star_ratio_upstream, T_star_ratio_downstream, rho_star_ratio_upstream, rho_star_ratio_downstream,\
    P0_star_ratio_upstream, P0_star_ratio_downstream, T0_star_ratio_upstream, T0_star_ratio_downstream, T0_downstream, T0_upstream, P0_star, T0_star")

variables = [R, gamma, q, M_upstream, M_downstream, T_upstream, T_downstream, P_upstream, P_downstream, rho_upstream, rho_downstream,\
                P0_upstream, P0_downstream, rho0_upstream, rho0_downstream, T_star, P_star, rho_star, a_upstream, a_downstream, a0_upstream,\
                a0_downstream, a_star]

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
    #P_upstream from P_star
    Eq(P_star * P_star_ratio_upstream, P_upstream),
    #P_downstream from P_star
    Eq(P_star * P_star_ratio_downstream, P_downstream),
    #T_upstream from T_star
    Eq(T_star * T_star_ratio_upstream, T_upstream),
    #T_downstream from T_star
    Eq(T_star * T_star_ratio_downstream, T_downstream),
    #rho_upstream from rho_star
    Eq(rho_star * rho_star_ratio_upstream, rho_upstream),
    #P_downstream from P_star
    Eq(rho_star * rho_star_ratio_downstream, rho_downstream),
    #P0_upstream from P0_star
    Eq(P0_star * P0_star_ratio_upstream, P0_upstream),
    #P0_downstream from P0_star
    Eq(P0_star * P0_star_ratio_downstream, P0_downstream),
    #T0_upstream from T0_star
    Eq(T0_star * T0_star_ratio_upstream, T0_upstream),
    #T0_downstream from P0_star
    Eq(T0_star * T0_star_ratio_downstream, T0_downstream),
    #T0_upstream from T_upstream
    Eq(T_upstream * (1 + (gamma - 1) / 2 * M_upstream ** 2), T0_upstream),
    #T0_downstream from T_downstream
    Eq(T_downstream * (1 + (gamma - 1) / 2 * M_downstream ** 2), T0_downstream),
    #P0_upstream from P_upstream
    Eq(P_upstream * (1 + (gamma - 1) / 2 * M_upstream ** 2) ** (gamma / (gamma - 1)), P0_upstream),
    #P0_downstream from P_downstream
    Eq(P_downstream * (1 + (gamma - 1) / 2 * M_downstream ** 2) ** (gamma / (gamma - 1)), P0_downstream),
    #rho0_upstream from rho_upstream
    Eq(rho_upstream * (1 + (gamma - 1) / 2 * M_upstream ** 2) ** (1 / (gamma - 1)), rho0_upstream),
    #rho0_downstream from rho_downstream
    Eq(rho_downstream * (1 + (gamma - 1) / 2 * M_downstream ** 2) ** (1 / (gamma - 1)), rho0_downstream),
    #T_star from a_star
    Eq(T0_downstream * (a_star / a0_downstream) ** 2, T_star),
    #T_star from a_star
    Eq(T0_upstream * (a_star / a0_upstream) ** 2, T_star),
    #T_star from gamma
    Eq(T0_downstream * 2 / (gamma + 1), T_star),
    #T_star from gamma
    Eq(T0_upstream * 2 / (gamma + 1), T_star),
    #a_star from gamma
    Eq(a0_downstream * sqrt(2 / (gamma + 1)), a_star),
    #a_star from gamma
    Eq(a0_upstream * sqrt(2 / (gamma + 1)), a_star),
    #P_star from P0_upstream
    Eq(P0_upstream * (2 / (gamma + 1)) ** (gamma / (gamma - 1)), P_star),
    #P_star from P0_downstream
    Eq(P0_downstream * (2 / (gamma + 1)) ** (gamma / (gamma - 1)), P_star),
    #rho_star from rho0_upstream
    Eq(rho0_upstream * (2 / (gamma + 1)) ** (1 / (gamma - 1)), rho_star),
    #rho_star from rho0_downstream
    Eq(rho0_downstream * (2 / (gamma + 1)) ** (1 / (gamma - 1)), rho_star),
]   

def adiabatic_relations(eqs, known_values):

    """
    Computes adiabatic flow properties given various inputs.
    Computes flow properties for upstream and downstream.
    """
    """ all_variables = set().union(*[eq.free_symbols for eq in eqs])
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
            break  # Stop if no new unknowns were solved """
    
    solution = solve(eqs, known_values)
    print(solution)

    return solution

def main():

    knowns = {}

    for var in variables:
        inputs = input(f"Enter {var} or press enter to skip: ")
        if inputs.strip():
            knowns[var] = float(inputs)
    
    print(knowns)
    
    #results = adiabatic_relations(equations, knowns)
    results = solve(equations, knowns[var])

    print(results)

    """ if results:
        print("\nSolved Values:")
        for var, value in results.items():
            print(f"{var} = {value:.4f}") """

        
if __name__ == "__main__":
    main()
