struct LogConfig
    occupation_number::Bool
    orbitals_energy::Bool
    stopping_criteria::Bool
    energy::Bool
    density::Bool

    function LogConfig(;occupation_number = false, orbitals_energy = false, stopping_criteria = true, energy = true, density = false)
        new(occupation_number, orbitals_energy, stopping_criteria, energy, density)
    end
end

struct LogBook
    config::LogConfig
    occupation_number_log
    orbitals_energy_log
    stopping_criteria_log
    energy_log
    density_log

    function LogBook(config, T)
        occupation_number_log = []                  # To see
        orbitals_energy_log   = []                  # To see
        stopping_criteria_log = T[]                 # To see
        energy_log            = T[]                 # To see
        density_log           = []                  # To see
        new(config, occupation_number_log, orbitals_energy_log, stopping_criteria_log, energy_log, density_log)
    end
end