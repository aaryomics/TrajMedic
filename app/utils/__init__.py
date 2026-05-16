# TrajMedic — utils package
from .helpers  import (load_universe, extract_metadata, get_time_array,
                        sample_frames, get_protein_selection,
                        get_calpha_selection, detect_jumps,
                        score_from_penalty, health_label)
from .plotting import (plot_rmsd, plot_rg, plot_com_drift,
                        plot_health_gauge,
                        plotly_rmsd, plotly_rg, plotly_com_drift)
from .report   import generate_txt_report, generate_html_report
