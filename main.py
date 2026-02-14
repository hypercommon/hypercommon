if __name__ == "__main__":
    '''from experiments.ring_lattice.thresholds_to_excel import run_experiment

    configs = [
        (40, 2),
        (40, 4),
        (40, 6),
        (40, 8),
        (40, 10),
        (40, 12),
        (40, 14),
        (40, 16),
        (40, 18)
    ]

    run_experiment(configs)'''

    from experiments.ring_lattice.rewiring_vs_truth import run_experiment
    configs = [(500, 16, 2)]

    run_experiment(
        configs=configs,
        avg_times=5,
        p_step=0.01
    )