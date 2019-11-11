
import numpy as np


class CallConsensusIntegrationSites:

    def __init__(self, region_size: int = 1000, minimum_observations: int = 5):
        self.region_size = region_size
        self.minimum_observations = minimum_observations

    def call_integration_sites(self, int_group: list = None):
        observed_integrations = self.get_region_calls(int_group)
        peak_regions = self.parse_peak_sites(observed_integrations)
        final_peaks = {}
        for int_site, supporting_calls in peak_regions.items():
            if len(supporting_calls[0]) > self.minimum_observations:
                final_peaks[int_site] = (int(np.mean(supporting_calls[0])), len(supporting_calls[0]),
                                         int(np.mean(supporting_calls[1])))
        return final_peaks

    def parse_peak_sites(self, observed_integrations):
        peak_regions = {}
        for chrom, sites in observed_integrations.items():
            sites.sort(key=lambda x: x[0])
            previous_sites, vector_sites = [], []
            for integration_site, vector_pos in sites:
                if not previous_sites:
                    previous_sites = [integration_site]
                    vector_sites = [vector_pos]
                else:
                    if integration_site - previous_sites[-1] > self.region_size:
                        peak_regions[f'{chrom}:{previous_sites[0]}-{previous_sites[-1]}'] = [previous_sites,
                                                                                             vector_sites]
                        previous_sites = [integration_site]
                        vector_sites = [vector_pos]
                    else:
                        previous_sites.append(integration_site)
                        vector_sites.append(vector_pos)
            if previous_sites:
                peak_regions[f'{chrom}:{previous_sites[0]}-{previous_sites[-1]}'] = [previous_sites, vector_sites]
        return peak_regions

    @staticmethod
    def get_region_calls(int_group):
        observed_integrations = {}
        for int_call in int_group:
            if int_call[2] not in observed_integrations:
                observed_integrations[int_call[2]] = [(int_call[3], int_call[-1])]
            else:
                observed_integrations[int_call[2]].append((int_call[3], int_call[-1]))
        return observed_integrations
