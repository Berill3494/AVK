#include <chrono>
#include <cmath>
#include <list>
#include <ncurses.h>
#include <random>
#include <ranges>
#include <thread>

using namespace std;
using namespace chrono_literals;

struct Random {
	random_device dev;
	mt19937 rng;
	uniform_int_distribution<> dist;
	Random() : rng(this->dev()), dist(0, 10000) {}

	double get(size_t m) {return this->dist(rng) % m;}
};

struct Vector {
	double x, y;

	void add(Vector const& other) {
		this->x += other.x;
		this->y += other.y;
	}
    void sub(Vector const& other) {
        this->x -= other.x;
        this->y -= other.y;
    }
    void scale(double scalar) {
        this->x *= scalar;
        this->y *= scalar;
    }
	void average(Vector const& other){
		add(other);
		this->x /= 2;
		this->y /= 2;
	}
    void length(double nl) {
        auto l = length();
        scale(nl/l);
    }
    void max(double l) {
        if(length() > l)
            length(l);
    }
    double distance(Vector const& other){
        return Vector(this->x-other.x, this->y-other.y).length(); 
    }
    double angle() const {
        return std::atan2(x, y);
    }
    double length() const {
        return std::sqrt(x*x + y*y);
    }
    Vector normalize() const {
        const auto l = length();
        return Vector(x/l, y/l);
    }
};

// ---------------------------------------------------------------------- //

namespace Constants {
    const double G = 6.67E-11; // gravitational constant
    const double c = 2.99792458; // speed of light (scaled)
    const double time_factor = .1; // how fast particles are travelling
    const double blackhole_mass = 500000000000;//E11;
    const Vector particle_velocity = {0, 0}; // initial velocity
    const size_t history_size = 10;
}

struct Particle {
	Vector	p;
	Vector  v;
    double  mass = 1;
	bool	show = true;
    std::list<std::pair<std::size_t, std::size_t>> history;

    void move() {
        std::pair<std::size_t, std::size_t> pp = {this->p.x, this->p.y};
        if(this->history.empty() || 
            this->history.front() != pp)
                this->history.push_front(pp);
        if(this->history.size() > Constants::history_size)
            this->history.resize(Constants::history_size);
        // move the particle (add velocity to position)
        this->p.add(this->v);     
    }

	void draw(WINDOW *w, Vector bh){
        if(this->p.x < 0 || this->p.y < 0)
            return;
 		if(!this->show)
			return;
        if(this->p.distance(bh) < 5)
            return;
        static string_view berill = "BERILL";
        std::size_t i = 1;
        for(auto const& p: history) {
            mvwaddch(w, p.second, p.first, i < berill.length() ? berill[i] : '.');
            ++i;
        }
        wattron(w, A_BOLD);
		mvwaddch(w, this->p.y, this->p.x, berill[0]);
        wattroff(w, A_BOLD);
	}
};


struct Blackhole {
	Vector p;
    double mass;
    const double event_horizon = 2*Constants::G*mass/(Constants::c*Constants::c);

    void draw(WINDOW *w) {
        if(this->p.x < 0 || this->p.y < 0)
            return;
        mvwaddch(w, this->p.y, this->p.x, ' ');
    }

	void pull(Particle &p){
        // calculate distance vector
		Vector dv = this->p;
        dv.sub(p.p);
		double distance = dv.length();
        double gravity = Constants::G*(this->mass*p.mass)/(distance*distance);

        // set to gravity
        dv.length(gravity);
        p.v.add(dv);
        // travel with the speed of light
        p.v.length(Constants::c * Constants::time_factor);

       if(distance < this->event_horizon)
			p.show = false;
	}
};


struct Screen {
	WINDOW *w = 0;
	size_t x, y = 0;

	void init() {
		initscr();
		curs_set(0);
        start_color();
        init_pair(1, COLOR_BLACK, COLOR_GREEN);
		getmaxyx(stdscr, this->y, this->x);
		this->w = newwin(y+1, x, 0, 0);
	}

	void refresh() {
		wrefresh(this->w);
	}
};




int main() {
   	Screen s;
	s.init();
	Random r;

	Blackhole b{{double(s.x/8), double(s.y/2)}, Constants::blackhole_mass};

    std::vector<Particle> ps;
    ps.reserve(s.x*s.y);
    mvaddstr(0, 0, string("screen size: " + std::to_string(s.x) + ":" + std::to_string(s.y)).c_str());
    s.refresh();


    const size_t height = s.y/2;
    for(auto i = 0u; i < height; ++i) {
        ps.push_back({{double(s.x-r.get(10)), double(i)}, Constants::particle_velocity});
    }


	for(;;) {
		wclear(s.w);
		for(auto &p: ps) {
            if(!p.show)
                continue;
			b.pull(p);
            p.move();
			p.draw(s.w, b.p);
		}
        b.draw(s.w);
		s.refresh();
		std::this_thread::sleep_for(50ms);
	}
	endwin();
}
