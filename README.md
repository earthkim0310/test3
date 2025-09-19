# Radial Velocity Exoplanet Simulator

A Streamlit-based interactive simulator for visualizing exoplanet radial velocity (RV) curves with real-time animation.

## Features

- 🌌 **Interactive RV Simulation**: Real-time visualization of exoplanet radial velocity curves
- ▶️ **Animation Controls**: Play/pause, speed control, and reset functionality
- ⚙️ **Parameter Adjustment**: 
  - Stellar mass (Ms)
  - Planet mass (Mp)
  - Orbital period (P)
  - Inclination angle (i)
  - Eccentricity (e)
  - Argument of periastron (ω)
- 📊 **Three Visualization Panels**:
  - RV curve with current time marker
  - Doppler shift of absorption lines
  - Orbital projection with current position

## Requirements

- streamlit>=1.28.0
- matplotlib>=3.7.0
- numpy>=1.24.0
- pandas>=2.0.0
- pillow>=9.0.0

## Local Installation

1. Clone the repository
2. Install dependencies:
   ```bash
   pip install -r requirements.txt
   ```
3. Run the application:
   ```bash
   streamlit run test3.py
   ```

## Streamlit Cloud Deployment

1. Push your code to GitHub
2. Go to [Streamlit Cloud](https://share.streamlit.io/)
3. Connect your GitHub repository
4. Select the main file: `test3.py`
5. Deploy!

## Usage

1. Adjust parameters in the sidebar
2. Use the play/pause button to start/stop animation
3. Adjust playback speed (0.25x to 4x)
4. Reset to return to initial state
5. Use the phase slider to set starting position

## Technical Details

- Uses matplotlib with Agg backend for cloud compatibility
- Implements Kepler's equation solver for elliptical orbits
- Real-time animation using Streamlit's rerun functionality
- Responsive design with three-column layout
