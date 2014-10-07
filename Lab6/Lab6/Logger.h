#pragma once
#include <fstream>
#include <functional>
namespace logger
{
    class log
    {
        std::ofstream m_fileStream;

    public:
        explicit log(char *string) : m_fileStream(std::ofstream(string))
        {
        }

        operator std::ostream & () { return m_fileStream; }

        #ifdef _DEBUG
        template <typename T>
        log &operator<<(const T &obj)
        {
            m_fileStream << obj;
            return *this;
        }
        log &operator <<(std::ostream& f(std::ostream&))
        {
            f(m_fileStream);
            return *this;
        }
        #else
        template <typename T>
        log &operator<<(const T &obj)
        {
            return *this;
        }
        #endif
    };
}