# function to find the closest time in a time sorted array of event times
# Times are MJD times, but any floating point array can be searched this way.
# HISTORY
# 25Aug10 GIL INitial version based on AI search

# initial version used date time structures
#from datetime import datetime
from bisect import bisect_left

def findTime(time_array, target_time):
    """
    Finds the closest time in a sorted array of irregularly sampled times
    to a given target time.

    Args:
$        time_array (list): A list of datetime objects, sorted in ascending order.
        time_array (list): List of MJD date+times values, sorted in ascending order.
        target_time (MJD f;pat): The target value object to find the closest match for.

    Returns:
        value closest to input value, and the index to the closest value
    """
    # check input quality
    if not time_array:
        return None, 0

    # Find the insertion point for the target_time
    # This gives the index where target_time would be inserted to maintain sort order.
    # If target_time is present, it returns the index of the first occurrence.
    idx = bisect_left(time_array, target_time)

    # Handle edge cases
    if idx == 0:  # Target time is before or equal to the first element
        return time_array[0], 0
    elif idx == len(time_array):  # Target time is after the last element
        return time_array[-1], len(time_array) - 1
    else:
        # Compare the element at idx and the element before idx
        # These are the two candidates for the closest time
        before_candidate = time_array[idx - 1]
        after_candidate = time_array[idx]

        diff_before = abs((target_time - before_candidate))
        diff_after = abs((target_time - after_candidate))

        if diff_after < diff_before:
            return after_candidate, idx
        else:
            return before_candidate, idx - 1

# Example Usage:
if __name__ == "__main__":
    # Create a sorted array of irregular sampled times
    times = [
        60800.1234, 60800.2345, 60800.3456, 60800.4567, 60800.5678, 60800.6789
    ]

    target = 60800.5
    closest, idx = findTime(times, target)
    print(f"The closest time to {target} is: {closest}, index is: {idx}")

    target_early = 60800.0
    closest_early, idx = findTime(times, target_early)
    print(f"The closest time to {target_early} is: {closest_early}, index is: {idx}")

    target_late = 60900.0
    closest_late, idx = findTime(times, target_late)
    print(f"The closest time to {target_late} is: {closest_late}, index is: {idx}")

    target_exact = 60800.2345
    closest_exact, idx = findTime(times, target_exact)
    print(f"The closest time to {target_exact} is: {closest_exact},index is: {idx}")

    empty_times = []
    closest_empty, idx = findTime(empty_times, target_exact)
    print(f"The closest time in an empty array is: {closest_empty}, index is: {idx}")
