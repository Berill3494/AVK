#include <chrono>
#include <cmath>
#include <iostream>
#include <list>
#include <ncurses.h>
#include <numbers>
#include <random>
#include <ranges>
#include <thread>
#include <tuple>

using namespace std;
using namespace chrono_literals;

struct Random {
	random_device dev;
	mt19937 rng;
	uniform_int_distribution<> dist;
	Random() : rng(this->dev()), dist(0, 10000) {}

	auto get(size_t m) {return this->dist(rng) % m;}
};

struct ScreenPosition {
    int x;
    int y;

    auto operator<=> (ScreenPosition const&) const = default;
};

struct Window {
    using ncurses_window = unique_ptr<WINDOW, decltype([](auto w){delwin(w);})>;
    ncurses_window w;
    Window(size_t x, size_t y, size_t px = 0, size_t py = 0) : w(newwin(y, x, py, px)) {
        nodelay(this->w.get(), 1);
    }

    auto get() const {
        return this->w.get();
    }

    auto refresh() {
		wrefresh(get());
	}

    auto get_max() const {
        struct {int xmax, ymax;} r;
        getmaxyx(get(), r.ymax, r.xmax);
        return r;
    }
    
    auto is_in(ScreenPosition const& p) const {
        if(p.x < 0 || p.y < 0)
            return false;
        auto [xmax, ymax] = get_max();
        if(p.x >= xmax || p.y >= ymax)
            return false;
        return true;
    }

    auto put(int x, int y, char ch) {
        mvwaddch(get(), y, x, ch);
    }
    auto put(int x, int y, string s) {
        mvwaddstr(get(), y, x, s.c_str());
    }

    auto bold(bool b) {
        b ? wattron(get(), A_BOLD) : wattroff(get(), A_BOLD);
    }
    
    auto clear() {
    	wclear(get());
    }

    auto key_pressed() const {
        int c = wgetch(get());
        return (c != ERR);
    }

    auto wait_for_key() {
        nodelay(get(), 0);
        int c = wgetch(get());
        nodelay(get(), 1);
        return c;
    }
};

struct Screen {
	int xmax, ymax = 0;
    Screen() {
        initscr();
        curs_set(0);
	    getmaxyx(stdscr, this->ymax, this->xmax);
    }
    ~Screen() {
        endwin();
    }

	auto get_window(size_t x, size_t y, size_t px = 0, size_t py = 0) {
        return Window(x, y, px, py);
	}
};

// ---------------------------------------------------------------------- //

struct Vector {
	double x, y;

    // adds another vector
	auto add(Vector const& other) {
		this->x += other.x;
		this->y += other.y;
	}
    // subtracts another vector
    auto sub(Vector const& other) {
        this->x -= other.x;
        this->y -= other.y;
    }
    // scales the vector with a scalar
    auto scale(double scalar) {
        this->x *= scalar;
        this->y *= scalar;
    }
    // scales the vector to the passed in magnitude
    auto magnitude(double nl) {
        auto l = magnitude();
        scale(nl/l);
    }
    // limits the magnitude of the vector
    auto limit(double l) {
        if(magnitude() > l)
            magnitude(l);
    }
    // creates a vector based on an angle in radians
    auto angle(double a) {
        this->y = 1;
        this->x = 1/tan(a);
    }
    // calculates the distance of two vectors
    auto distance(Vector const& other) {
        return Vector(this->x-other.x, this->y-other.y).magnitude(); 
    }
    // returns the angle in radians
    auto angle() const -> double {
        return atan2(y, x);
    }
    // returns the magnitude (length) of a vector
    auto magnitude() const -> double {
        // Pythaghoras
        return sqrt(x*x + y*y);
    }
    // returns the normalized version of this vector (with magnitude 1)
    auto normalize() const {
        const auto l = magnitude();
        return Vector(x/l, y/l);
    }

    auto to_screen(double scale) -> ScreenPosition const {
        return {int(this->x*scale), int(this->y*scale)};
    }

    static auto from_screen(ScreenPosition const& sp, double scale) -> Vector {
        // make the calculation in double
        double x = sp.x, y = sp.y;
        return {x/scale, y/scale};
    }
};

// calculates the difference of two vectors
auto difference(Vector const& v0, Vector const& v1) {
    Vector v = v0;
    v.sub(v1);
    return v;
}

// ---------------------------------------------------------------------- //

namespace Constants {
    // gravitational constant (https://en.wikipedia.org/wiki/Gravitational_constant)
    const double G = 6.6743E-11;
    // m/s speed of light in vacuum
    const double c = 299792458;

    // M (solar mass) this determines the size and force of the black hole Sagittarius A* (https://en.wikipedia.org/wiki/Sagittarius_A*)
    const double blackhole_mass = 4.297E6;
    const double solar_mass = 1.989E30;

    // this determines the frame rate of the movement
    const double time_factor = 1;
    // scaling of the universe
    const double scale = 1E-9;
    // initial direction of the particles
    const Vector particle_velocity = {-1, 0};
    // how long the trail is for the particles
    const size_t history_size = 10;
}

struct Particle {
    // the position of the particle in our geometry (https://en.wikipedia.org/wiki/2D_computer_graphics)
	Vector	p;
    // speed of the particle (velocity)
	Vector  v;
    // particles have zero mass, but for simplicity we will work with 1
    double  mass = 1;
    // stores information if the particle needs to be calculated or is already eaten by the blackhole
	bool	show = true;
    // stores the old screen positions of the particle to be able to draw the "trail"
    // it stores the position in integers not in double as its easier to check for equality and at the
    // end the screen positions are integers
    list<ScreenPosition> history;


    Vector rv;
    double last_r;
    double last_f;
    Vector last_v;

    auto move(Vector const& f) {
        ScreenPosition last = p.to_screen(Constants::scale);
        // only store the last position if it is different from the last one stored
        // this ensures that even if the particle ramained in the same position, it will have a trail
        if(this->history.empty() || this->history.front() != last)
                this->history.push_front(last);
        // limit the history to the trail size
        if(this->history.size() > Constants::history_size)
            this->history.resize(Constants::history_size);

        last_v = this->v;
        this->v.add(f);
        // ensure that the particle always moves with the speed of light
        this->v.magnitude(Constants::c);
        Vector pv = this->v;
        // scale it to our time frame
        pv.scale(Constants::time_factor);
        // move the particle (add velocity to position)
        rv = pv;
        this->p.add(pv);     
    }

	auto draw(Window &w) {
        auto sp = this->p.to_screen(Constants::scale);
        w.put(0,  9, "vlast: " + to_string(last_v.x) + ":" + to_string(last_v.y) + ", mag: " + to_string(last_v.magnitude()));
        w.put(0, 10, "v: " + to_string(rv.x) + ":" + to_string(rv.y) + ", mag: " + to_string(rv.magnitude()));
        w.put(0, 11, "p: " + to_string(p.x) + ":" + to_string(p.y) + ", distance: " + to_string(last_r));
        w.put(0, 12, "d: " + to_string(last_r));
        w.put(0, 13, "f: " + to_string(last_f));
        w.put(0, 14, "p: " + to_string(sp.x) + ":" + to_string(sp.y) + ", show: " + to_string(show));

        if(!this->show || !w.is_in(sp))
            return;
        static string_view berill = "BERILL";
        // we start from 1 as the first character will always be printed even when there is no tail
        std::size_t i = 1;
        for(auto const& p: history) {
            w.put(p.x, p.y, i < berill.length() ? berill[i] : '.');
            ++i;
        }
        w.put(0, 1, "angle: " + to_string(this->v.angle()) + ", v.angle: " + to_string(this->v.angle()));
        w.bold(true);
		w.put(sp.x, sp.y, berill[0]);
        w.bold(false);
	}
};


struct Blackhole {
	Vector p;
    double mass; // in solar mass
    // calculate the Schwarzschild radius (event horizon) of the blackhole (https://en.wikipedia.org/wiki/Schwarzschild_radius)
    // r = 2 × G × M/c*c
    const double r = 2*Constants::G*mass*Constants::solar_mass/(Constants::c*Constants::c);

    auto draw(Window &w) {
        auto sp = this->p.to_screen(Constants::scale);
        if(!w.is_in(sp))
            return;

        const int R = this->r*Constants::scale;
        // calculate the bounding square for the circle, add 1 to make sure the whole circle is included
        // as there might be rounding mistakes
        const int N = 2*R+1;

        // exit if we have set a weird large universe
        auto [xmax, ymax] = w.get_max();
        if(N>min(xmax, ymax)) {
            cout << "Blackhole too big\n";
            exit(-1);
        }

        // simple filled circle raster algorithm, which goes through each
        // pixed in the bounding box of the circle and fills it, it its inside
        // the circle (distance is smaller than radius)
        for(int i = 0; i < N; i++) {
            for(int j = 0; j < N; j++) {
                int x = i-R;
                int y = j-R;

                // Pythagoras to see if the point is inside the triangle to the rim of the circle
                if (x*x + y*y <= R*R+1)
                    w.put(sp.x-N/2+i, sp.y-N/2+j, 'o');
            }
        }
    }
    auto newtonian_pull(Particle &p) {
        // calculate the distance vector between the two center points
        // the magnitude
        Vector towards_blackhole = difference(this->p, p.p);
        // the distance from the particle to the blackhole is the length of the vector
        // Newton's law of universal gravity (https://en.wikipedia.org/wiki/Newton%27s_law_of_universal_gravitation)
        //    F = G*m0*m1/(r*r)
        // where F the gravitaion force, G is the gravitation constant, m0-m1 are the masses of the two objects, r is the distance between m0 and m1
        auto r = towards_blackhole.magnitude();
        auto F = Constants::G*this->mass*Constants::solar_mass*p.mass/(r*r);
        p.last_f = F;
        p.last_r = r*r;
        // we already have the direction of the force from the vector, now we just need to limit its force to F by setting its length
        towards_blackhole.magnitude(F);
        return towards_blackhole;
    }

	auto proper_pull(Particle &p){
        /*// calculate distance vector
		Vector dv = this->p;
        dv.sub(p.p);
        auto angle = dv.angle();
		double distance = dv.magnitude();
        double gravity = Constants::G*(this->mass*p.mass)/(distance*distance);

        //double delta_angle = gravity*(Constants::time_factor/Constants::c) * std::sin(p.angle - angle);
        //delta_angle /= std::abs(1.-2.*Constants::G*this->mass/(distance*Constants::c*Constants::c));
        //p.angle += delta_angle;
        //p.angle += std::sin(p.angle - angle);

        p.v.angle(angle);
        p.v.magnitude(Constants::c);
        
        if(distance < this->r)
            p.show = false;*/
    }
};

// Declare it static, so we can exit() and still restore the terminal states
static Screen s;


int main() {
	Random r;
    Window w = s.get_window(s.xmax, s.ymax);
    auto [xm, ym] = w.get_max();

	Blackhole b{Vector::from_screen({s.xmax/8, s.ymax/2}, Constants::scale), Constants::blackhole_mass};

    vector<Particle> ps;

    const int height = s.ymax/2;
    for(int i = 0; i < height; ++i) {
        ps.push_back({Vector::from_screen({s.xmax-int(r.get(10)), i}, Constants::scale), Constants::particle_velocity});
    }

	for(;;) {
        w.clear();
        w.put(0, 0, "resolution: " + to_string(xm) + ":" + to_string(ym));
		for(auto &p: ps) {
            // do not even calculate dead particles
            if(!p.show)
                continue;
			auto f = b.newtonian_pull(p);
            p.move(f);
			p.draw(w);
		}
        b.draw(w);
		w.refresh();
        if(w.key_pressed())
            w.wait_for_key();
		this_thread::sleep_for(100ms);
	}
}
