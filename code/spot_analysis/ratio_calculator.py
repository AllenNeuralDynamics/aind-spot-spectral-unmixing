import torch
import numpy as np
from pathlib import Path
from typing import Tuple
import os
from .config import Config

class RatioCalculator:
    def __init__(self):
        os.environ['CUDA_LAUNCH_BLOCKING'] = '1'  # Enable CUDA launch blocking
        os.environ['TORCH_USE_CUDA_DSA'] = '1'

        self.config = Config
        # self.device = torch.device('cuda:0' if torch.cuda.is_available() else 'cpu')#
        self.channels = self.config.get_round_channels().keys()
        
    # def objective_fn(self, r: torch.Tensor, subset: torch.Tensor, L1: float) -> torch.Tensor:
    #     """Calculate objective function for ratio optimization"""
    #     r = r / torch.norm(r, dim=0)
    #     n_cam = len(self.config.get_round_channels())
        
    #     dot_products = torch.tile(subset @ r, (n_cam, 1, 1))
    #     ys = torch.tile(torch.unsqueeze(r, 1), (1, subset.shape[0], 1))
    #     xs = torch.tile(
    #         torch.unsqueeze(torch.transpose(subset, 0, 1), 2),
    #         (1, 1, n_cam)
    #     )
        
    #     return (
    #         torch.sum(torch.min(torch.norm(dot_products * ys - xs, dim=0), dim=1)[0]) +
    #         L1 * torch.sum(torch.abs(r))
    #     )
    
    # def calculate_ratios(self, intensity_data: np.ndarray, ratio_path: Path) -> np.ndarray:
    #     """Calculate or load channel ratios"""
    #     os.environ['CUDA_LAUNCH_BLOCKING'] = '1'  # Enable CUDA launch blocking
    #     os.environ['TORCH_USE_CUDA_DSA'] = '1'
    #     torch.backends.cuda.matmul.allow_tf32 = False  # Disable TensorFloat-32 (TF32) for better error reporting

    #     if ratio_path.exists():
    #         return np.loadtxt(ratio_path).T
        
    #     # Initialize parameters
    #     n_cam = len(self.config.get_round_channels())
    #     initial = np.eye(len(self.channels))
    #     n_sub = int(self.config.N_SUBSET * self.config.FRAC_SAMPLED)
        
    #     # Setup tensors
    #     r_gpu = torch.from_numpy(initial).cuda(0).requires_grad_()
    #     # data_gpu = torch.from_numpy(intensity_data).to(self.device).double()
    #     data_gpu = torch.from_numpy(np.array(intensity_data)).cuda(0).double()
        
    #     # Optimization loop
    #     loss_hist = torch.zeros(self.config.EPOCHS)
    #     r_hist = torch.zeros((self.config.EPOCHS, n_cam, n_cam))
        
    #     for i in range(self.config.EPOCHS):
    #         if not i % self.config.RESAMPLE_ITER:
    #             sub_gpu = data_gpu[np.random.choice(
    #                 self.config.N_SUBSET,
    #                 n_sub,
    #                 replace=False
    #             )]
            
    #         loss = self.objective_fn(r_gpu.double(), sub_gpu, self.config.L1)
    #         loss.backward()
            
    #         loss_hist[i] = self.objective_fn(
    #             r_gpu.clone().double(),
    #             data_gpu,
    #             self.config.L1
    #         ).data
    #         r_hist[i] = r_gpu.clone().double()
            
    #         r_gpu.data -= self.config.LEARNING_RATE * r_gpu.grad.data
    #         r_gpu.data = torch.div(r_gpu.data, torch.norm(r_gpu.data))
    #         r_gpu.grad = None
        
    #     # Get optimized ratios
    #     optimized = r_hist[np.argmin(loss_hist)].detach().numpy()
    #     optimized = optimized / np.linalg.norm(optimized, axis=0)
        
    #     # Save ratios
    #     np.savetxt(
    #         ratio_path,
    #         100 * optimized.T / optimized.max(0)[..., None],
    #         delimiter='\t',
    #         fmt='%d'
    #     )
        
    #     return optimized

    def subset_spots_df(self, thresh_spots):
        if len(thresh_spots)> self.config.N_SUBSET:
            thresh_spots_subsetted = thresh_spots[::int(len(thresh_spots)/self.config.N_SUBSET)]

        self.config.N_SUBSET = len(thresh_spots_subsetted)
        return thresh_spots_subsetted

    def calculate_ratios(self, intensity_data: np.ndarray, ratio_location: Path) -> np.ndarray:
        if not os.path.exists(ratio_location):
            # Max normalizing the ratios
            def norm(ratios): 
                return (100 * ratios / ratios.max()).astype(int)

            # Linear least squares objective function
            def objective_fn(r, subset, L1): 
                n_cam = len(self.config.get_round_channels())
                r = r / torch.norm(r, dim=0)
                dot_products = torch.tile(subset @ r, (n_cam, 1, 1))
                ys = torch.tile(torch.unsqueeze(r,1), (1, subset.shape[0], 1))
                xs = torch.tile(torch.unsqueeze(torch.transpose(subset, 0, 1),2), (1, 1, n_cam))
                return torch.sum(torch.min(torch.norm(dot_products * ys - xs, dim=0), dim=1)[0]) + self.config.L1 * torch.sum(torch.abs(r)) 
            n_cam = len(self.config.get_round_channels())
            initial = np.eye((n_cam)) # Initial ratios matrices for each channel

            n_sub = np.int32(self.config.N_SUBSET * self.config.FRAC_SAMPLED)
            os.environ['CUDA_LAUNCH_BLOCKING'] = '1'  # Enable CUDA launch blocking
            os.environ['TORCH_USE_CUDA_DSA'] = '1'
            torch.backends.cuda.matmul.allow_tf32 = False  # Disable TensorFloat-32 (TF32) for better error reporting
            # initial = np.random.rand(n_cam, n_dye) # for checking that initialized performance is much better than with random lines
            r_gpu = torch.from_numpy(initial).cuda(0).requires_grad_()
            thresh_spots_df_subsetted = self.subset_spots_df(intensity_data)
            data_gpu = torch.from_numpy(np.array(thresh_spots_df_subsetted)).cuda(0).double()

            loss_hist = torch.zeros(self.config.EPOCHS)
            r_hist = torch.zeros((self.config.EPOCHS, n_cam, n_cam))

            for i in range(self.config.EPOCHS):
                if not i%self.config.RESAMPLE_ITER:
                    sub_gpu = data_gpu[np.random.choice(self.config.N_SUBSET, n_sub, replace=False)]
                if not i % 1000 and i or i == 1: print(i, loss_hist[i-1], flush=True)
                loss = objective_fn(r_gpu.double(), sub_gpu, self.config.L1)
                loss.backward()
                loss_hist[i] = objective_fn(r_gpu.clone().double(), data_gpu, self.config.L1).data
                # loss_hist[i] = objective_fn(r_gpu.clone().double(), sub_gpu, self.config.L1).data

                r_hist[i] = r_gpu.clone().double()
                r_gpu.data -= self.config.LEARNING_RATE * r_gpu.grad.data
                r_gpu.data = torch.div(r_gpu.data, torch.norm(r_gpu.data))
                r_gpu.grad = None

            optimized = r_hist[np.argmin(loss_hist)].detach().numpy()
            optimized = optimized/np.linalg.norm(optimized, axis=0)
            np.savetxt(ratio_location, 100 * optimized.T / optimized.max(0)[..., None], delimiter='\t', fmt='%d')
            print('original')
            for x in initial.T: print(norm(x))
            print('optimized')
            for x in optimized.T: print(norm(x))
            print('error', loss_hist.min()/len(intensity_data))
        ratios = np.loadtxt(ratio_location).T
        return ratios