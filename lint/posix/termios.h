
/*  termios.h  --  POSIX header
    Copyright 1992, Gimpel Software, All Rights Reserved
 */

#if defined(_lint) && !defined(__TERMIOS_H)
#define __TERMIOS_H

typedef  unsigned  tcflag_t;
typedef  unsigned  cc_t;
typedef  unsigned  speed_t;

struct termios {
	tcflag_t   c_iflag;
	tcflag_t   c_oflag;
	tcflag_t   c_cflag;
	tcflag_t   c_lflag;
	cc_t       c_cc;
};

speed_t cfgetispeed(const struct termios * );
speed_t cfgetospeed(const struct termios * );
int cfsetispeed(struct termios *, speed_t );
int cfsetospeed(struct termios *, speed_t );
int tcdrain(int );
int tcflow(int, int );
int tcflush(int, int );
int tcgetattr(int, struct termios * );
int tcsetattr(int, int, const struct termios * );
int tcsendbreak(int, int );

#define B0          1
#define B110        2
#define B1200       3
#define B134        4
#define B150        5
#define B1800       6
#define B19200      7
#define B200        8
#define B2400       9
#define B300       10
#define B38400     11
#define B4800      12
#define B50        13
#define B600       14
#define B75        15
#define B9600      16
#define BRKINT     17
#define CLOCAL     18
#define CREAD      19
#define CS5        20
#define CS6        21
#define CS7        22
#define CS8        23
#define CSIZE      24
#define CSTOPB     25
#define ECHO       26
#define ECHOE      27
#define ECHOK      28
#define ECHONL     29
#define HUPCL      30
#define ICANON     31
#define ICRNL      32
#define IEXTEN     33
#define IGNBRK     34
#define IGNCR      35
#define IGNPAR     36
#define INLCR      37
#define INPCK      38
#define ISIG       39
#define ISTRIP     40
#define IXOFF      41
#define IXON       42
#define NCCS       43
#define NOFLSH     44
#define OPOST      45
#define PARENB     46
#define PARMRK     47
#define PARODD     48
#define TCIFLUSH   49
#define TCIOFF     50
#define TCIOFLUSH  51
#define TCION      52
#define TCOFLUSH   53
#define TCOOFF     54
#define TCOON      55
#define TCSADRAIN  56
#define TCSAFLUSH  57
#define TCSANOW    58
#define TOSTOP     59
#define VEOF       60
#define VEOL       61
#define VERASE     62
#define VINTR      63
#define VKILL      64
#define VMIN       65
#define VQUIT      66
#define VSTART     67
#define VSTOP      68
#define VSUSP      69
#define VTIME      70

#endif
