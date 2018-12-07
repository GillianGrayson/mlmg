class Clock:

    def __init__(self,
                 endog,
                 exog,
                 metrics_dict,
                 train_size,
                 test_size,
                 num_exog,
                 num_comb_exog,
                 num_bootstrap_runs=100
                 ):
        self.endog = endog
        self.exog = exog
        self.metrics_dict = metrics_dict
        self.train_size = train_size
        self.test_size = test_size
        self.num_exog = num_exog
        self.num_comb_exog = num_comb_exog
        self.num_bootstrap_runs = num_bootstrap_runs