#include <chrono>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <string>
#include <thread>

void show_progress_bar(std::ostream& os, int time,
                       std::string message, char symbol = '*')
{
    static const auto bar_length = 70;
    // not including the percentage figure and spaces

    if (message.length() >= bar_length) {
        os << message << '\n';
        message.clear();
    } else {
        message += " ";
    }

    const auto progress_level = 100.0 / (bar_length - message.length());

    std::cout << message;

    for (double percentage = 0; percentage <= 100; percentage += progress_level) {
        message += symbol;
        os << "\r [" << std::setw(3) << static_cast<int>(percentage) << "%] "
           << message << std::flush;
        std::this_thread::sleep_for(std::chrono::milliseconds(time));
    }
    os << "\n\n";
}


int main()
{
    show_progress_bar(std::clog, 100, "progress", '#');
}