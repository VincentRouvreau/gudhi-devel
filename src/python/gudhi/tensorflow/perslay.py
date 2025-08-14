import torch
import torch.nn as nn
import torch.nn.functional as F
import math

class GridPerslayWeight(nn.Module):
    def __init__(self, grid, grid_bnds):
        super().__init__()
        self.grid = nn.Parameter(torch.tensor(grid, dtype=torch.float32))
        self.grid_bnds = torch.tensor(grid_bnds, dtype=torch.float32)

    def forward(self, diagrams):
        # diagrams: list of [num_pts, 2] tensors, batch size n
        grid_shape = self.grid.shape
        min_x, max_x = self.grid_bnds[0]
        min_y, max_y = self.grid_bnds[1]
        # For each diagram in batch
        weights = []
        for d in diagrams:
            # d: [num_pts, 2]
            ids_x = ((grid_shape[0] * (d[:, 0] - min_x) / (max_x - min_x))).long().clamp(0, grid_shape[0]-1)
            ids_y = ((grid_shape[1] * (d[:, 1] - min_y) / (max_y - min_y))).long().clamp(0, grid_shape[1]-1)
            w = self.grid[ids_x, ids_y]
            weights.append(w)
        return weights

class GaussianMixturePerslayWeight(nn.Module):
    def __init__(self, gaussians):
        super().__init__()
        self.W = nn.Parameter(torch.tensor(gaussians, dtype=torch.float32))  # shape [4, n]

    def forward(self, diagrams):
        # diagrams: list of [num_pts, 2] tensors
        means = self.W[:2, :].transpose(0,1)  # [n, 2]
        variances = self.W[2:, :].transpose(0,1)  # [n, 2]
        weights = []
        for d in diagrams:
            # d: [num_pts, 2]
            d_exp = d.unsqueeze(1)  # [num_pts, 1, 2]
            m_exp = means.unsqueeze(0)  # [1, n, 2]
            v_exp = variances.unsqueeze(0)  # [1, n, 2]
            dists = ((d_exp - m_exp) ** 2) / (v_exp ** 2)
            w = torch.exp(-dists.sum(dim=2)).sum(dim=1)
            weights.append(w)
        return weights

class PowerPerslayWeight(nn.Module):
    def __init__(self, constant, power):
        super().__init__()
        self.constant = nn.Parameter(torch.tensor(constant, dtype=torch.float32))
        self.power = power

    def forward(self, diagrams):
        weights = []
        for d in diagrams:
            dist = torch.abs(d[:, 1] - d[:, 0])
            w = self.constant * torch.pow(dist, self.power)
            weights.append(w)
        return weights

class GaussianPerslayPhi(nn.Module):
    def __init__(self, image_size, image_bnds, variance):
        super().__init__()
        self.image_size = image_size
        self.image_bnds = image_bnds
        self.variance = nn.Parameter(torch.tensor(variance, dtype=torch.float32))

    def forward(self, diagrams):
        # diagrams: list of [num_pts, 2] tensors
        step = [(self.image_bnds[i][1] - self.image_bnds[i][0]) / self.image_size[i] for i in range(2)]
        coords = [torch.arange(self.image_bnds[i][0], self.image_bnds[i][1], step[i]) for i in range(2)]
        M = torch.meshgrid(coords[0], coords[1], indexing='ij')
        mu = torch.stack([tens for tens in M], dim=0)  # [2, n_x, n_y]
        output_list = []
        output_shape = M[0].shape + (1,)
        for d in diagrams:
            t = torch.as_tensor(d, dtype=torch.float32)  # Convert numpy to tensor
            d_d = torch.stack([t[:, 0], t[:, 1] - t[:, 0]], dim=1)  # [num_pts, 2]
            # Broadcast for gaussian eval
            for _ in range(2):
                d_d = d_d.unsqueeze(-1)
            mu_exp = mu.unsqueeze(0)  # [1, 2, n_x, n_y]
            dists = ((d_d - mu_exp) ** 2) / (2 * self.variance ** 2)
            gauss = torch.exp(-dists.sum(dim=1)) / (2 * math.pi * (self.variance ** 2))
            output = gauss.unsqueeze(-1)
            output_list.append(output)
        return output_list, output_shape

class TentPerslayPhi(nn.Module):
    def __init__(self, samples):
        super().__init__()
        self.samples = nn.Parameter(torch.tensor(samples, dtype=torch.float32))

    def forward(self, diagrams):
        output_list = []
        output_shape = self.samples.shape
        for d in diagrams:
            xs = d[:, 0:1]  # [num_pts,1]
            ys = d[:, 1:2]
            samples_d = self.samples.unsqueeze(0).unsqueeze(0)  # [1,1,num_samples]
            val = 0.5 * (ys - xs) - torch.abs(samples_d - 0.5 * (ys + xs))
            output = torch.maximum(val, torch.tensor(0.0))
            output_list.append(output.squeeze(1))  # [num_pts, num_samples]
        return output_list, output_shape

class FlatPerslayPhi(nn.Module):
    def __init__(self, samples, theta):
        super().__init__()
        self.samples = nn.Parameter(torch.tensor(samples, dtype=torch.float32))
        self.theta = nn.Parameter(torch.tensor(theta, dtype=torch.float32))

    def forward(self, diagrams):
        output_list = []
        output_shape = self.samples.shape
        for d in diagrams:
            xs = d[:, 0:1]
            ys = d[:, 1:2]
            samples_d = self.samples.unsqueeze(0).unsqueeze(0)
            val = 0.5 * (ys - xs) - torch.abs(samples_d - 0.5 * (ys + xs))
            output = 1.0 / (1.0 + torch.exp(-self.theta * val))
            output_list.append(output.squeeze(1))
        return output_list, output_shape

class Perslay(nn.Module):
    def __init__(self, weight, phi, perm_op, rho):
        super().__init__()
        self.weight = weight
        self.phi = phi
        self.perm_op = perm_op
        self.rho = rho

    def forward(self, diagrams):
        # diagrams: list of [num_pts, 2] tensors
        vector_list, dim = self.phi(diagrams)
        weight_list = self.weight(diagrams)
        output_list = []
        for v, w in zip(vector_list, weight_list):
            # expand weight to match v's shape
            for _ in range(len(dim)):
                w = w.unsqueeze(-1)
            vw = v * w
            if isinstance(self.perm_op, str) and self.perm_op[:3] == "top":
                k = int(self.perm_op[3:])
                # pad vw to fixed length
                vw_pad = F.pad(vw, (0,0,0, max(0, k*dim[0]-vw.shape[0])), value=-1e10)
                vw_reshaped = vw_pad.transpose(0,1).reshape(-1)
                topk, _ = torch.topk(vw_reshaped, k*dim[0])
                out = topk
            else:
                out = self.perm_op(vw, dim=0)
            out2 = self.rho(out)
            output_list.append(out2)
        # Stack outputs for batch
        return torch.stack(output_list)
