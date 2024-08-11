
void calculateJac(Eigen::Matrix<double, 3, 4>& Jac){
    double roll, pitch, yaw, q1, q2, q3,q4;
    double dx1, dy1, dz1, dx2, dy2, dz2 , dx3, dy3, dz3, dx4, dy4, dz4, dx5, dy5, dz5;

    Jac << -(dy2*cos(pitch)*cos(q1)*cos(yaw)) - dy3*cos(pitch)*cos(q1)*cos(q2)*cos(yaw) - dy4*cos(pitch)*cos(q1)*cos(q2)*cos(yaw) - dy5*cos(pitch)*cos(q1)*cos(q2)*cos(yaw) + dz3*cos(pitch)*cos(q1)*cos(yaw)*sin(q2) +
           dx2*(-(cos(pitch)*cos(yaw)*sin(q1)) + cos(q1)*(cos(yaw)*sin(pitch)*sin(roll) - cos(roll)*sin(yaw))) + dx3*(-(cos(pitch)*cos(yaw)*sin(q1)) + cos(q1)*(cos(yaw)*sin(pitch)*sin(roll) - cos(roll)*sin(yaw))) +
           dx4*(-(cos(pitch)*cos(q1)*cos(yaw)*sin(q2)*sin(q3)) + cos(q3)*(-(cos(pitch)*cos(yaw)*sin(q1)) + cos(q1)*(cos(yaw)*sin(pitch)*sin(roll) - cos(roll)*sin(yaw)))) +
           dz4*(cos(pitch)*cos(q1)*cos(q3)*cos(yaw)*sin(q2) + sin(q3)*(-(cos(pitch)*cos(yaw)*sin(q1)) + cos(q1)*(cos(yaw)*sin(pitch)*sin(roll) - cos(roll)*sin(yaw)))) +
           dz5*(-(sin(q3 - q4)*(-(cos(pitch)*cos(q1)*cos(yaw)*sin(q2)*sin(q3)) + cos(q3)*(-(cos(pitch)*cos(yaw)*sin(q1)) + cos(q1)*(cos(yaw)*sin(pitch)*sin(roll) - cos(roll)*sin(yaw))))) +
                cos(q3 - q4)*(cos(pitch)*cos(q1)*cos(q3)*cos(yaw)*sin(q2) + sin(q3)*(-(cos(pitch)*cos(yaw)*sin(q1)) + cos(q1)*(cos(yaw)*sin(pitch)*sin(roll) - cos(roll)*sin(yaw))))) +
           dx5*(cos(q3 - q4)*(-(cos(pitch)*cos(q1)*cos(yaw)*sin(q2)*sin(q3)) + cos(q3)*(-(cos(pitch)*cos(yaw)*sin(q1)) + cos(q1)*(cos(yaw)*sin(pitch)*sin(roll) - cos(roll)*sin(yaw)))) +
                sin(q3 - q4)*(cos(pitch)*cos(q1)*cos(q3)*cos(yaw)*sin(q2) + sin(q3)*(-(cos(pitch)*cos(yaw)*sin(q1)) + cos(q1)*(cos(yaw)*sin(pitch)*sin(roll) - cos(roll)*sin(yaw))))),
            -(dy2*sin(q2)*(cos(yaw)*sin(pitch)*sin(roll) - cos(roll)*sin(yaw))) + dz3*(pow(sin(q2),2)*(cos(yaw)*sin(pitch)*sin(roll) - cos(roll)*sin(yaw)) - sin(q2)*(cos(roll)*cos(yaw)*sin(pitch) + sin(roll)*sin(yaw)) -
                                                                                       cos(q2)*(-(cos(pitch)*cos(yaw)*sin(q1)) + cos(q2)*(cos(yaw)*sin(pitch)*sin(roll) - cos(roll)*sin(yaw)))) +
            dz4*cos(q3)*(pow(sin(q2),2)*(cos(yaw)*sin(pitch)*sin(roll) - cos(roll)*sin(yaw)) - sin(q2)*(cos(roll)*cos(yaw)*sin(pitch) + sin(roll)*sin(yaw)) -
                         cos(q2)*(-(cos(pitch)*cos(yaw)*sin(q1)) + cos(q2)*(cos(yaw)*sin(pitch)*sin(roll) - cos(roll)*sin(yaw)))) -
            dx4*sin(q3)*(pow(sin(q2),2)*(cos(yaw)*sin(pitch)*sin(roll) - cos(roll)*sin(yaw)) - sin(q2)*(cos(roll)*cos(yaw)*sin(pitch) + sin(roll)*sin(yaw)) -
                         cos(q2)*(-(cos(pitch)*cos(yaw)*sin(q1)) + cos(q2)*(cos(yaw)*sin(pitch)*sin(roll) - cos(roll)*sin(yaw)))) +
            dy3*(-(cos(q2)*sin(q2)*(cos(yaw)*sin(pitch)*sin(roll) - cos(roll)*sin(yaw))) + cos(q2)*(cos(roll)*cos(yaw)*sin(pitch) + sin(roll)*sin(yaw)) - sin(q2)*(-(cos(pitch)*cos(yaw)*sin(q1)) + cos(q2)*(cos(yaw)*sin(pitch)*sin(roll) - cos(roll)*sin(yaw)))) +
            dy4*(-(cos(q2)*sin(q2)*(cos(yaw)*sin(pitch)*sin(roll) - cos(roll)*sin(yaw))) + cos(q2)*(cos(roll)*cos(yaw)*sin(pitch) + sin(roll)*sin(yaw)) - sin(q2)*(-(cos(pitch)*cos(yaw)*sin(q1)) + cos(q2)*(cos(yaw)*sin(pitch)*sin(roll) - cos(roll)*sin(yaw)))) +
            dy5*(-(cos(q2)*sin(q2)*(cos(yaw)*sin(pitch)*sin(roll) - cos(roll)*sin(yaw))) + cos(q2)*(cos(roll)*cos(yaw)*sin(pitch) + sin(roll)*sin(yaw)) - sin(q2)*(-(cos(pitch)*cos(yaw)*sin(q1)) + cos(q2)*(cos(yaw)*sin(pitch)*sin(roll) - cos(roll)*sin(yaw)))) +
            dx5*(-(cos(q3 - q4)*sin(q3)*(pow(sin(q2),2)*(cos(yaw)*sin(pitch)*sin(roll) - cos(roll)*sin(yaw)) - sin(q2)*(cos(roll)*cos(yaw)*sin(pitch) + sin(roll)*sin(yaw)) -
                                         cos(q2)*(-(cos(pitch)*cos(yaw)*sin(q1)) + cos(q2)*(cos(yaw)*sin(pitch)*sin(roll) - cos(roll)*sin(yaw))))) +
                 cos(q3)*sin(q3 - q4)*(pow(sin(q2),2)*(cos(yaw)*sin(pitch)*sin(roll) - cos(roll)*sin(yaw)) - sin(q2)*(cos(roll)*cos(yaw)*sin(pitch) + sin(roll)*sin(yaw)) -
                                       cos(q2)*(-(cos(pitch)*cos(yaw)*sin(q1)) + cos(q2)*(cos(yaw)*sin(pitch)*sin(roll) - cos(roll)*sin(yaw))))) +
            dz5*(cos(q3)*cos(q3 - q4)*(pow(sin(q2),2)*(cos(yaw)*sin(pitch)*sin(roll) - cos(roll)*sin(yaw)) - sin(q2)*(cos(roll)*cos(yaw)*sin(pitch) + sin(roll)*sin(yaw)) -
                                       cos(q2)*(-(cos(pitch)*cos(yaw)*sin(q1)) + cos(q2)*(cos(yaw)*sin(pitch)*sin(roll) - cos(roll)*sin(yaw)))) +
                 sin(q3)*sin(q3 - q4)*(pow(sin(q2),2)*(cos(yaw)*sin(pitch)*sin(roll) - cos(roll)*sin(yaw)) - sin(q2)*(cos(roll)*cos(yaw)*sin(pitch) + sin(roll)*sin(yaw)) -
                                       cos(q2)*(-(cos(pitch)*cos(yaw)*sin(q1)) + cos(q2)*(cos(yaw)*sin(pitch)*sin(roll) - cos(roll)*sin(yaw))))),
            dx4*(-(sin(q3)*(cos(pitch)*cos(q1)*cos(yaw) + sin(q1)*(cos(yaw)*sin(pitch)*sin(roll) - cos(roll)*sin(yaw)))) -
                 cos(q3)*(cos(q2)*(cos(roll)*cos(yaw)*sin(pitch) + sin(roll)*sin(yaw)) - sin(q2)*(-(cos(pitch)*cos(yaw)*sin(q1)) + cos(q2)*(cos(yaw)*sin(pitch)*sin(roll) - cos(roll)*sin(yaw))))) +
            dz4*(cos(q3)*(cos(pitch)*cos(q1)*cos(yaw) + sin(q1)*(cos(yaw)*sin(pitch)*sin(roll) - cos(roll)*sin(yaw))) -
                 sin(q3)*(cos(q2)*(cos(roll)*cos(yaw)*sin(pitch) + sin(roll)*sin(yaw)) - sin(q2)*(-(cos(pitch)*cos(yaw)*sin(q1)) + cos(q2)*(cos(yaw)*sin(pitch)*sin(roll) - cos(roll)*sin(yaw))))) +
            dx5*(cos(q3 - q4)*(-(sin(q3)*(cos(pitch)*cos(q1)*cos(yaw) + sin(q1)*(cos(yaw)*sin(pitch)*sin(roll) - cos(roll)*sin(yaw)))) -
                               cos(q3)*(cos(q2)*(cos(roll)*cos(yaw)*sin(pitch) + sin(roll)*sin(yaw)) - sin(q2)*(-(cos(pitch)*cos(yaw)*sin(q1)) + cos(q2)*(cos(yaw)*sin(pitch)*sin(roll) - cos(roll)*sin(yaw))))) +
                 cos(q3 - q4)*(sin(q3)*(cos(pitch)*cos(q1)*cos(yaw) + sin(q1)*(cos(yaw)*sin(pitch)*sin(roll) - cos(roll)*sin(yaw))) +
                               cos(q3)*(cos(q2)*(cos(roll)*cos(yaw)*sin(pitch) + sin(roll)*sin(yaw)) - sin(q2)*(-(cos(pitch)*cos(yaw)*sin(q1)) + cos(q2)*(cos(yaw)*sin(pitch)*sin(roll) - cos(roll)*sin(yaw)))))) +
            dz5*(-(sin(q3 - q4)*(-(sin(q3)*(cos(pitch)*cos(q1)*cos(yaw) + sin(q1)*(cos(yaw)*sin(pitch)*sin(roll) - cos(roll)*sin(yaw)))) -
                                 cos(q3)*(cos(q2)*(cos(roll)*cos(yaw)*sin(pitch) + sin(roll)*sin(yaw)) - sin(q2)*(-(cos(pitch)*cos(yaw)*sin(q1)) + cos(q2)*(cos(yaw)*sin(pitch)*sin(roll) - cos(roll)*sin(yaw)))))) -
                 sin(q3 - q4)*(sin(q3)*(cos(pitch)*cos(q1)*cos(yaw) + sin(q1)*(cos(yaw)*sin(pitch)*sin(roll) - cos(roll)*sin(yaw))) +
                               cos(q3)*(cos(q2)*(cos(roll)*cos(yaw)*sin(pitch) + sin(roll)*sin(yaw)) - sin(q2)*(-(cos(pitch)*cos(yaw)*sin(q1)) + cos(q2)*(cos(yaw)*sin(pitch)*sin(roll) - cos(roll)*sin(yaw)))))),
            dz5*(sin(q3 - q4)*(sin(q3)*(cos(pitch)*cos(q1)*cos(yaw) + sin(q1)*(cos(yaw)*sin(pitch)*sin(roll) - cos(roll)*sin(yaw))) +
                               cos(q3)*(cos(q2)*(cos(roll)*cos(yaw)*sin(pitch) + sin(roll)*sin(yaw)) - sin(q2)*(-(cos(pitch)*cos(yaw)*sin(q1)) + cos(q2)*(cos(yaw)*sin(pitch)*sin(roll) - cos(roll)*sin(yaw))))) +
                 cos(q3 - q4)*(cos(q3)*(cos(pitch)*cos(q1)*cos(yaw) + sin(q1)*(cos(yaw)*sin(pitch)*sin(roll) - cos(roll)*sin(yaw))) -
                               sin(q3)*(cos(q2)*(cos(roll)*cos(yaw)*sin(pitch) + sin(roll)*sin(yaw)) - sin(q2)*(-(cos(pitch)*cos(yaw)*sin(q1)) + cos(q2)*(cos(yaw)*sin(pitch)*sin(roll) - cos(roll)*sin(yaw)))))) +
            dx5*(-(cos(q3 - q4)*(sin(q3)*(cos(pitch)*cos(q1)*cos(yaw) + sin(q1)*(cos(yaw)*sin(pitch)*sin(roll) - cos(roll)*sin(yaw))) +
                                 cos(q3)*(cos(q2)*(cos(roll)*cos(yaw)*sin(pitch) + sin(roll)*sin(yaw)) - sin(q2)*(-(cos(pitch)*cos(yaw)*sin(q1)) + cos(q2)*(cos(yaw)*sin(pitch)*sin(roll) - cos(roll)*sin(yaw)))))) +
                 sin(q3 - q4)*(cos(q3)*(cos(pitch)*cos(q1)*cos(yaw) + sin(q1)*(cos(yaw)*sin(pitch)*sin(roll) - cos(roll)*sin(yaw))) -
                               sin(q3)*(cos(q2)*(cos(roll)*cos(yaw)*sin(pitch) + sin(roll)*sin(yaw)) - sin(q2)*(-(cos(pitch)*cos(yaw)*sin(q1)) + cos(q2)*(cos(yaw)*sin(pitch)*sin(roll) - cos(roll)*sin(yaw)))))),
            -(dy2*cos(pitch)*cos(q1)*sin(yaw)) - dy3*cos(pitch)*cos(q1)*cos(q2)*sin(yaw) - dy4*cos(pitch)*cos(q1)*cos(q2)*sin(yaw) - dy5*cos(pitch)*cos(q1)*cos(q2)*sin(yaw) + dz3*cos(pitch)*cos(q1)*sin(q2)*sin(yaw) +
            dx2*(-(cos(pitch)*sin(q1)*sin(yaw)) + cos(q1)*(cos(roll)*cos(yaw) + sin(pitch)*sin(roll)*sin(yaw))) + dx3*(-(cos(pitch)*sin(q1)*sin(yaw)) + cos(q1)*(cos(roll)*cos(yaw) + sin(pitch)*sin(roll)*sin(yaw))) +
            dx4*(-(cos(pitch)*cos(q1)*sin(q2)*sin(q3)*sin(yaw)) + cos(q3)*(-(cos(pitch)*sin(q1)*sin(yaw)) + cos(q1)*(cos(roll)*cos(yaw) + sin(pitch)*sin(roll)*sin(yaw)))) +
            dz4*(cos(pitch)*cos(q1)*cos(q3)*sin(q2)*sin(yaw) + sin(q3)*(-(cos(pitch)*sin(q1)*sin(yaw)) + cos(q1)*(cos(roll)*cos(yaw) + sin(pitch)*sin(roll)*sin(yaw)))) +
            dz5*(-(sin(q3 - q4)*(-(cos(pitch)*cos(q1)*sin(q2)*sin(q3)*sin(yaw)) + cos(q3)*(-(cos(pitch)*sin(q1)*sin(yaw)) + cos(q1)*(cos(roll)*cos(yaw) + sin(pitch)*sin(roll)*sin(yaw))))) +
                 cos(q3 - q4)*(cos(pitch)*cos(q1)*cos(q3)*sin(q2)*sin(yaw) + sin(q3)*(-(cos(pitch)*sin(q1)*sin(yaw)) + cos(q1)*(cos(roll)*cos(yaw) + sin(pitch)*sin(roll)*sin(yaw))))) +
            dx5*(cos(q3 - q4)*(-(cos(pitch)*cos(q1)*sin(q2)*sin(q3)*sin(yaw)) + cos(q3)*(-(cos(pitch)*sin(q1)*sin(yaw)) + cos(q1)*(cos(roll)*cos(yaw) + sin(pitch)*sin(roll)*sin(yaw)))) +
                 sin(q3 - q4)*(cos(pitch)*cos(q1)*cos(q3)*sin(q2)*sin(yaw) + sin(q3)*(-(cos(pitch)*sin(q1)*sin(yaw)) + cos(q1)*(cos(roll)*cos(yaw) + sin(pitch)*sin(roll)*sin(yaw))))),
            -(dy2*sin(q2)*(cos(roll)*cos(yaw) + sin(pitch)*sin(roll)*sin(yaw))) + dz3*(-(sin(q2)*(-(cos(yaw)*sin(roll)) + cos(roll)*sin(pitch)*sin(yaw))) + pow(sin(q2),2)*(cos(roll)*cos(yaw) + sin(pitch)*sin(roll)*sin(yaw)) -
                                                                                       cos(q2)*(-(cos(pitch)*sin(q1)*sin(yaw)) + cos(q2)*(cos(roll)*cos(yaw) + sin(pitch)*sin(roll)*sin(yaw)))) +
            dz4*cos(q3)*(-(sin(q2)*(-(cos(yaw)*sin(roll)) + cos(roll)*sin(pitch)*sin(yaw))) + pow(sin(q2),2)*(cos(roll)*cos(yaw) + sin(pitch)*sin(roll)*sin(yaw)) -
                         cos(q2)*(-(cos(pitch)*sin(q1)*sin(yaw)) + cos(q2)*(cos(roll)*cos(yaw) + sin(pitch)*sin(roll)*sin(yaw)))) -
            dx4*sin(q3)*(-(sin(q2)*(-(cos(yaw)*sin(roll)) + cos(roll)*sin(pitch)*sin(yaw))) + pow(sin(q2),2)*(cos(roll)*cos(yaw) + sin(pitch)*sin(roll)*sin(yaw)) -
                         cos(q2)*(-(cos(pitch)*sin(q1)*sin(yaw)) + cos(q2)*(cos(roll)*cos(yaw) + sin(pitch)*sin(roll)*sin(yaw)))) +
            dy3*(cos(q2)*(-(cos(yaw)*sin(roll)) + cos(roll)*sin(pitch)*sin(yaw)) - cos(q2)*sin(q2)*(cos(roll)*cos(yaw) + sin(pitch)*sin(roll)*sin(yaw)) - sin(q2)*(-(cos(pitch)*sin(q1)*sin(yaw)) + cos(q2)*(cos(roll)*cos(yaw) + sin(pitch)*sin(roll)*sin(yaw)))) +
            dy4*(cos(q2)*(-(cos(yaw)*sin(roll)) + cos(roll)*sin(pitch)*sin(yaw)) - cos(q2)*sin(q2)*(cos(roll)*cos(yaw) + sin(pitch)*sin(roll)*sin(yaw)) - sin(q2)*(-(cos(pitch)*sin(q1)*sin(yaw)) + cos(q2)*(cos(roll)*cos(yaw) + sin(pitch)*sin(roll)*sin(yaw)))) +
            dy5*(cos(q2)*(-(cos(yaw)*sin(roll)) + cos(roll)*sin(pitch)*sin(yaw)) - cos(q2)*sin(q2)*(cos(roll)*cos(yaw) + sin(pitch)*sin(roll)*sin(yaw)) - sin(q2)*(-(cos(pitch)*sin(q1)*sin(yaw)) + cos(q2)*(cos(roll)*cos(yaw) + sin(pitch)*sin(roll)*sin(yaw)))) +
            dx5*(-(cos(q3 - q4)*sin(q3)*(-(sin(q2)*(-(cos(yaw)*sin(roll)) + cos(roll)*sin(pitch)*sin(yaw))) + pow(sin(q2),2)*(cos(roll)*cos(yaw) + sin(pitch)*sin(roll)*sin(yaw)) -
                                         cos(q2)*(-(cos(pitch)*sin(q1)*sin(yaw)) + cos(q2)*(cos(roll)*cos(yaw) + sin(pitch)*sin(roll)*sin(yaw))))) +
                 cos(q3)*sin(q3 - q4)*(-(sin(q2)*(-(cos(yaw)*sin(roll)) + cos(roll)*sin(pitch)*sin(yaw))) + pow(sin(q2),2)*(cos(roll)*cos(yaw) + sin(pitch)*sin(roll)*sin(yaw)) -
                                       cos(q2)*(-(cos(pitch)*sin(q1)*sin(yaw)) + cos(q2)*(cos(roll)*cos(yaw) + sin(pitch)*sin(roll)*sin(yaw))))) +
            dz5*(cos(q3)*cos(q3 - q4)*(-(sin(q2)*(-(cos(yaw)*sin(roll)) + cos(roll)*sin(pitch)*sin(yaw))) + pow(sin(q2),2)*(cos(roll)*cos(yaw) + sin(pitch)*sin(roll)*sin(yaw)) -
                                       cos(q2)*(-(cos(pitch)*sin(q1)*sin(yaw)) + cos(q2)*(cos(roll)*cos(yaw) + sin(pitch)*sin(roll)*sin(yaw)))) +
                 sin(q3)*sin(q3 - q4)*(-(sin(q2)*(-(cos(yaw)*sin(roll)) + cos(roll)*sin(pitch)*sin(yaw))) + pow(sin(q2),2)*(cos(roll)*cos(yaw) + sin(pitch)*sin(roll)*sin(yaw)) -
                                       cos(q2)*(-(cos(pitch)*sin(q1)*sin(yaw)) + cos(q2)*(cos(roll)*cos(yaw) + sin(pitch)*sin(roll)*sin(yaw))))),
            dx4*(-(sin(q3)*(cos(pitch)*cos(q1)*sin(yaw) + sin(q1)*(cos(roll)*cos(yaw) + sin(pitch)*sin(roll)*sin(yaw)))) -
                 cos(q3)*(cos(q2)*(-(cos(yaw)*sin(roll)) + cos(roll)*sin(pitch)*sin(yaw)) - sin(q2)*(-(cos(pitch)*sin(q1)*sin(yaw)) + cos(q2)*(cos(roll)*cos(yaw) + sin(pitch)*sin(roll)*sin(yaw))))) +
            dz4*(cos(q3)*(cos(pitch)*cos(q1)*sin(yaw) + sin(q1)*(cos(roll)*cos(yaw) + sin(pitch)*sin(roll)*sin(yaw))) -
                 sin(q3)*(cos(q2)*(-(cos(yaw)*sin(roll)) + cos(roll)*sin(pitch)*sin(yaw)) - sin(q2)*(-(cos(pitch)*sin(q1)*sin(yaw)) + cos(q2)*(cos(roll)*cos(yaw) + sin(pitch)*sin(roll)*sin(yaw))))) +
            dx5*(cos(q3 - q4)*(-(sin(q3)*(cos(pitch)*cos(q1)*sin(yaw) + sin(q1)*(cos(roll)*cos(yaw) + sin(pitch)*sin(roll)*sin(yaw)))) -
                               cos(q3)*(cos(q2)*(-(cos(yaw)*sin(roll)) + cos(roll)*sin(pitch)*sin(yaw)) - sin(q2)*(-(cos(pitch)*sin(q1)*sin(yaw)) + cos(q2)*(cos(roll)*cos(yaw) + sin(pitch)*sin(roll)*sin(yaw))))) +
                 cos(q3 - q4)*(sin(q3)*(cos(pitch)*cos(q1)*sin(yaw) + sin(q1)*(cos(roll)*cos(yaw) + sin(pitch)*sin(roll)*sin(yaw))) +
                               cos(q3)*(cos(q2)*(-(cos(yaw)*sin(roll)) + cos(roll)*sin(pitch)*sin(yaw)) - sin(q2)*(-(cos(pitch)*sin(q1)*sin(yaw)) + cos(q2)*(cos(roll)*cos(yaw) + sin(pitch)*sin(roll)*sin(yaw)))))) +
            dz5*(-(sin(q3 - q4)*(-(sin(q3)*(cos(pitch)*cos(q1)*sin(yaw) + sin(q1)*(cos(roll)*cos(yaw) + sin(pitch)*sin(roll)*sin(yaw)))) -
                                 cos(q3)*(cos(q2)*(-(cos(yaw)*sin(roll)) + cos(roll)*sin(pitch)*sin(yaw)) - sin(q2)*(-(cos(pitch)*sin(q1)*sin(yaw)) + cos(q2)*(cos(roll)*cos(yaw) + sin(pitch)*sin(roll)*sin(yaw)))))) -
                 sin(q3 - q4)*(sin(q3)*(cos(pitch)*cos(q1)*sin(yaw) + sin(q1)*(cos(roll)*cos(yaw) + sin(pitch)*sin(roll)*sin(yaw))) +
                               cos(q3)*(cos(q2)*(-(cos(yaw)*sin(roll)) + cos(roll)*sin(pitch)*sin(yaw)) - sin(q2)*(-(cos(pitch)*sin(q1)*sin(yaw)) + cos(q2)*(cos(roll)*cos(yaw) + sin(pitch)*sin(roll)*sin(yaw)))))),
            dz5*(sin(q3 - q4)*(sin(q3)*(cos(pitch)*cos(q1)*sin(yaw) + sin(q1)*(cos(roll)*cos(yaw) + sin(pitch)*sin(roll)*sin(yaw))) +
                               cos(q3)*(cos(q2)*(-(cos(yaw)*sin(roll)) + cos(roll)*sin(pitch)*sin(yaw)) - sin(q2)*(-(cos(pitch)*sin(q1)*sin(yaw)) + cos(q2)*(cos(roll)*cos(yaw) + sin(pitch)*sin(roll)*sin(yaw))))) +
                 cos(q3 - q4)*(cos(q3)*(cos(pitch)*cos(q1)*sin(yaw) + sin(q1)*(cos(roll)*cos(yaw) + sin(pitch)*sin(roll)*sin(yaw))) -
                               sin(q3)*(cos(q2)*(-(cos(yaw)*sin(roll)) + cos(roll)*sin(pitch)*sin(yaw)) - sin(q2)*(-(cos(pitch)*sin(q1)*sin(yaw)) + cos(q2)*(cos(roll)*cos(yaw) + sin(pitch)*sin(roll)*sin(yaw)))))) +
            dx5*(-(cos(q3 - q4)*(sin(q3)*(cos(pitch)*cos(q1)*sin(yaw) + sin(q1)*(cos(roll)*cos(yaw) + sin(pitch)*sin(roll)*sin(yaw))) +
                                 cos(q3)*(cos(q2)*(-(cos(yaw)*sin(roll)) + cos(roll)*sin(pitch)*sin(yaw)) - sin(q2)*(-(cos(pitch)*sin(q1)*sin(yaw)) + cos(q2)*(cos(roll)*cos(yaw) + sin(pitch)*sin(roll)*sin(yaw)))))) +
                 sin(q3 - q4)*(cos(q3)*(cos(pitch)*cos(q1)*sin(yaw) + sin(q1)*(cos(roll)*cos(yaw) + sin(pitch)*sin(roll)*sin(yaw))) -
                               sin(q3)*(cos(q2)*(-(cos(yaw)*sin(roll)) + cos(roll)*sin(pitch)*sin(yaw)) - sin(q2)*(-(cos(pitch)*sin(q1)*sin(yaw)) + cos(q2)*(cos(roll)*cos(yaw) + sin(pitch)*sin(roll)*sin(yaw)))))),
            dy2*cos(q1)*sin(pitch) + dy3*cos(q1)*cos(q2)*sin(pitch) + dy4*cos(q1)*cos(q2)*sin(pitch) + dy5*cos(q1)*cos(q2)*sin(pitch) - dz3*cos(q1)*sin(pitch)*sin(q2) + dx2*(sin(pitch)*sin(q1) + cos(pitch)*cos(q1)*sin(roll)) +
            dx3*(sin(pitch)*sin(q1) + cos(pitch)*cos(q1)*sin(roll)) + dx4*(cos(q1)*sin(pitch)*sin(q2)*sin(q3) + cos(q3)*(sin(pitch)*sin(q1) + cos(pitch)*cos(q1)*sin(roll))) +
            dz4*(-(cos(q1)*cos(q3)*sin(pitch)*sin(q2)) + sin(q3)*(sin(pitch)*sin(q1) + cos(pitch)*cos(q1)*sin(roll))) +
            dz5*(-(sin(q3 - q4)*(cos(q1)*sin(pitch)*sin(q2)*sin(q3) + cos(q3)*(sin(pitch)*sin(q1) + cos(pitch)*cos(q1)*sin(roll)))) + cos(q3 - q4)*(-(cos(q1)*cos(q3)*sin(pitch)*sin(q2)) + sin(q3)*(sin(pitch)*sin(q1) + cos(pitch)*cos(q1)*sin(roll)))) +
            dx5*(cos(q3 - q4)*(cos(q1)*sin(pitch)*sin(q2)*sin(q3) + cos(q3)*(sin(pitch)*sin(q1) + cos(pitch)*cos(q1)*sin(roll))) + sin(q3 - q4)*(-(cos(q1)*cos(q3)*sin(pitch)*sin(q2)) + sin(q3)*(sin(pitch)*sin(q1) + cos(pitch)*cos(q1)*sin(roll)))),
            -(dy2*cos(pitch)*sin(q2)*sin(roll)) + dz3*(-(cos(pitch)*cos(roll)*sin(q2)) + cos(pitch)*pow(sin(q2),2)*sin(roll) - cos(q2)*(sin(pitch)*sin(q1) + cos(pitch)*cos(q2)*sin(roll))) +
            dz4*cos(q3)*(-(cos(pitch)*cos(roll)*sin(q2)) + cos(pitch)*pow(sin(q2),2)*sin(roll) - cos(q2)*(sin(pitch)*sin(q1) + cos(pitch)*cos(q2)*sin(roll))) -
            dx4*sin(q3)*(-(cos(pitch)*cos(roll)*sin(q2)) + cos(pitch)*pow(sin(q2),2)*sin(roll) - cos(q2)*(sin(pitch)*sin(q1) + cos(pitch)*cos(q2)*sin(roll))) +
            dy3*(cos(pitch)*cos(q2)*cos(roll) - cos(pitch)*cos(q2)*sin(q2)*sin(roll) - sin(q2)*(sin(pitch)*sin(q1) + cos(pitch)*cos(q2)*sin(roll))) +
            dy4*(cos(pitch)*cos(q2)*cos(roll) - cos(pitch)*cos(q2)*sin(q2)*sin(roll) - sin(q2)*(sin(pitch)*sin(q1) + cos(pitch)*cos(q2)*sin(roll))) +
            dy5*(cos(pitch)*cos(q2)*cos(roll) - cos(pitch)*cos(q2)*sin(q2)*sin(roll) - sin(q2)*(sin(pitch)*sin(q1) + cos(pitch)*cos(q2)*sin(roll))) +
            dx5*(-(cos(q3 - q4)*sin(q3)*(-(cos(pitch)*cos(roll)*sin(q2)) + cos(pitch)*pow(sin(q2),2)*sin(roll) - cos(q2)*(sin(pitch)*sin(q1) + cos(pitch)*cos(q2)*sin(roll)))) +
                 cos(q3)*sin(q3 - q4)*(-(cos(pitch)*cos(roll)*sin(q2)) + cos(pitch)*pow(sin(q2),2)*sin(roll) - cos(q2)*(sin(pitch)*sin(q1) + cos(pitch)*cos(q2)*sin(roll)))) +
            dz5*(cos(q3)*cos(q3 - q4)*(-(cos(pitch)*cos(roll)*sin(q2)) + cos(pitch)*pow(sin(q2),2)*sin(roll) - cos(q2)*(sin(pitch)*sin(q1) + cos(pitch)*cos(q2)*sin(roll))) +
                 sin(q3)*sin(q3 - q4)*(-(cos(pitch)*cos(roll)*sin(q2)) + cos(pitch)*pow(sin(q2),2)*sin(roll) - cos(q2)*(sin(pitch)*sin(q1) + cos(pitch)*cos(q2)*sin(roll)))),
            dx4*(-(sin(q3)*(-(cos(q1)*sin(pitch)) + cos(pitch)*sin(q1)*sin(roll))) - cos(q3)*(cos(pitch)*cos(q2)*cos(roll) - sin(q2)*(sin(pitch)*sin(q1) + cos(pitch)*cos(q2)*sin(roll)))) +
            dz4*(cos(q3)*(-(cos(q1)*sin(pitch)) + cos(pitch)*sin(q1)*sin(roll)) - sin(q3)*(cos(pitch)*cos(q2)*cos(roll) - sin(q2)*(sin(pitch)*sin(q1) + cos(pitch)*cos(q2)*sin(roll)))) +
            dx5*(cos(q3 - q4)*(-(sin(q3)*(-(cos(q1)*sin(pitch)) + cos(pitch)*sin(q1)*sin(roll))) - cos(q3)*(cos(pitch)*cos(q2)*cos(roll) - sin(q2)*(sin(pitch)*sin(q1) + cos(pitch)*cos(q2)*sin(roll)))) +
                 cos(q3 - q4)*(sin(q3)*(-(cos(q1)*sin(pitch)) + cos(pitch)*sin(q1)*sin(roll)) + cos(q3)*(cos(pitch)*cos(q2)*cos(roll) - sin(q2)*(sin(pitch)*sin(q1) + cos(pitch)*cos(q2)*sin(roll))))) +
            dz5*(-(sin(q3 - q4)*(-(sin(q3)*(-(cos(q1)*sin(pitch)) + cos(pitch)*sin(q1)*sin(roll))) - cos(q3)*(cos(pitch)*cos(q2)*cos(roll) - sin(q2)*(sin(pitch)*sin(q1) + cos(pitch)*cos(q2)*sin(roll))))) -
                 sin(q3 - q4)*(sin(q3)*(-(cos(q1)*sin(pitch)) + cos(pitch)*sin(q1)*sin(roll)) + cos(q3)*(cos(pitch)*cos(q2)*cos(roll) - sin(q2)*(sin(pitch)*sin(q1) + cos(pitch)*cos(q2)*sin(roll))))),
            dz5*(sin(q3 - q4)*(sin(q3)*(-(cos(q1)*sin(pitch)) + cos(pitch)*sin(q1)*sin(roll)) + cos(q3)*(cos(pitch)*cos(q2)*cos(roll) - sin(q2)*(sin(pitch)*sin(q1) + cos(pitch)*cos(q2)*sin(roll)))) +
                 cos(q3 - q4)*(cos(q3)*(-(cos(q1)*sin(pitch)) + cos(pitch)*sin(q1)*sin(roll)) - sin(q3)*(cos(pitch)*cos(q2)*cos(roll) - sin(q2)*(sin(pitch)*sin(q1) + cos(pitch)*cos(q2)*sin(roll))))) +
            dx5*(-(cos(q3 - q4)*(sin(q3)*(-(cos(q1)*sin(pitch)) + cos(pitch)*sin(q1)*sin(roll)) + cos(q3)*(cos(pitch)*cos(q2)*cos(roll) - sin(q2)*(sin(pitch)*sin(q1) + cos(pitch)*cos(q2)*sin(roll))))) +
                 sin(q3 - q4)*(cos(q3)*(-(cos(q1)*sin(pitch)) + cos(pitch)*sin(q1)*sin(roll)) - sin(q3)*(cos(pitch)*cos(q2)*cos(roll) - sin(q2)*(sin(pitch)*sin(q1) + cos(pitch)*cos(q2)*sin(roll)))));
}