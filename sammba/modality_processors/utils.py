import sys
import time
from nilearn.datasets.utils import _format_time


def _iterate_and_show_progress(iterator, total_number):
    t0 = time.time()
    iterator_elements = []
    initial_n = 0
    for n_so_far in range(total_number):
        iterator_elements.append(iterator.next())

        # Estimate remaining download time
        total_percent = float(n_so_far) / total_number
    
        current_download_size = n_so_far - initial_n
        n_remaining = total_number - n_so_far
        dt = time.time() - t0
        computation_rate = current_download_size / max(1e-8, float(dt))
        # Minimum rate of 0.01 bytes/s, to avoid dividing by zero.
        time_remaining = n_remaining / max(0.01, computation_rate)
    
        # Trailing whitespace is to erase extra char when message length
        # varies
        sys.stderr.write(
            "\rComputed %d of %d (%.1f%%, %s remaining)"
            % (n_so_far, total_number, total_percent * 100,
            _format_time(time_remaining)))
        initial_n += 1
    return iterator_elements