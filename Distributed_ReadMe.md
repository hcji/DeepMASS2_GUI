# **DeepMASS2 Distributed System**

## 1. Environment Dependencies(Central Server、Subserver、User Interface)

1. Clone the repository and navigate into it:

   ```bash
   git clone https://github.com/hcji/DeepMASS2_GUI.git
   cd DeepMASS2_GUI
   ```

2. Create and activate the Conda environment:

   ```bash
   conda env create -f environment.yml
   conda activate deepmass
   ```

3. Install additional dependencies:

   ```bash
   pip install cryptography
   ```

4. Model and database files (only required when running `client.py` on subservers)

Download the following files and place them in the `model/` directory: [dependent data](https://github.com/hcji/DeepMASS2_GUI/releases/tag/v0.99.1)

```bash
Ms2Vec_allGNPSnegative.hdf5
Ms2Vec_allGNPSnegative.hdf5.syn1neg.npy
Ms2Vec_allGNPSnegative.hdf5.wv.vectors.npy
Ms2Vec_allGNPSpositive.hdf5
Ms2Vec_allGNPSpositive.hdf5.syn1neg.npy
Ms2Vec_allGNPSpositive.hdf5.wv.vectors.npy
```

Reference Database (Optional)

- If you **do not have a custom database**, place **all** the following files in the `data/` directory: [dependent data](https://github.com/hcji/DeepMASS2_GUI/releases/tag/v0.99.1)

  ```bash
  DeepMassStructureDB-v1.1.csv
  references_index_negative_spec2vec.bin
  references_index_positive_spec2vec.bin
  references_spectrums_negative.pickle
  references_spectrums_positive.pickle
  ```

- If you **use a custom database**:

  - Keep `DeepMassStructureDB-v1.1.csv` unchanged in the `data/` directory.
  - Replace the `.pickle` and `.bin` files with your own data (you may keep the original filenames).
  - Or update the paths in `client.py` to point to your custom files.

> [!IMPORTANT]
>
> Refer to [`vectorize_reference_by_ms2vec.py`](https://github.com/hcji/DeepMASS2_Data_Processing/blob/master/Scripts/training_models/vectorize_reference_by_ms2vec.py) to regenerate each `.bin` file from its corresponding `.pickle` file.

***Note**: Step 4 is only required when running `client.py`.*

------

## 2. Deployment Architecture

| Role           | Script                    | Startup Example                                              |
| -------------- | ------------------------- | ------------------------------------------------------------ |
| Central Server | `server.py`               | `python server.py <SERVER_IP> <REG_PORT> <FILE_PORT>`        |
| Subserver      | `client.py`               | `python client.py <SERVER_IP> <REG_PORT> <FILE_PORT> <CLIENT_IP> <CLIENT_FILE_PORT>` |
| User Interface | `DistributedDeepMASS2.py` | In the GUI prompt, enter `<SERVER_IP>:<FILE_PORT>` and run   |

- **CENTRAL_SERVER_IP**: IP address of the central server (e.g., `192.168.1.10`)
- **SUBSERVER_IP**:   IP address of each subserver (e.g., `192.168.1.3`, `192.168.1.4`, …)
- **REG_PORT**:       Port for subserver registration (e.g., `5001`)
- **FILE_PORT**:      Port for user file uploads on the central server (e.g., `5002`)
- **CLIENT_FILE_PORT**: Port on each subserver for listening to file requests (e.g., `6001`, `6002`, …)

------

## 3. Usage Steps

### 1. Check Local IP

- **Linux**:

  ```bash
  ip addr show   # or `hostname -I`
  ```

- **Windows**:

  ```bash
  ipconfig
  ```

*Look for the `inet` / `IPv4 Address` entry to find your LAN IP.*

### 2. Start the Central Server

**Run on the central server (e.g., at IP `192.168.1.10`):**

```bash
cd DeepMASS2_GUI/server
python server.py 192.168.1.10 5001 5002
```

- Listening for registrations: `5001`
- Listening for user file uploads: `5002`

### 3. Start Subservers (Multiple Instances Supported)

**Run on the subserver (e.g., at IP `192.168.1.3`):**

```bash
cd DeepMASS2_GUI
python client.py 192.168.1.10 5001 5002 192.168.1.3 6001
```

- `192.168.1.10`: Central server IP
- `5001`: Registration port
- `5002`: User file upload port
- `192.168.1.3`: This subserver’s IP
- `6001`: File service listening port on this subserver

*Subservers will automatically register with the central server upon startup.*

### 4. Start the User Interface

```bash
cd DeepMASS2_GUI
python DistributedDeepMASS2.py
```

1. Click **IP Address**, then enter the central server’s IP and file-port in the format
    `<SERVER_IP>:<FILE_PORT>` (for example, `192.168.1.10:5002`).
2. Click **Open** to upload your mass spectrometry file.
3. Click **Run DeepMASS** to encrypt and upload your file to the central server, which then distributes it to all subservers.
4. To save results, click **Save** to export them as CSV files.
